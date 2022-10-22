#include <chrono>
#include <complex>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <utility>
#include <vector>

enum MPI_TAGS { MPI_TAG_END_PROCESS, MPI_TAG_MANDELBROT_TASK };

class Task {
public:
    virtual void run() = 0;
    virtual void sendTask(unsigned int to) = 0;
    virtual void getTask(unsigned int from) = 0;
    virtual void getResult(unsigned int from) = 0;
    virtual void sendResult(unsigned int to) = 0;

    void finalizeProcess(unsigned int to){
        MPI_Send(nullptr, 0, MPI_INT, to, MPI_TAG_END_PROCESS, MPI_COMM_WORLD);
    }
};

class MandelbrotTask : public Task {
private:
    std::complex<float> c;
    unsigned int N_iter;
    bool ans;
    std::pair<int,int> id;

public:
    void setC(std::complex<float> c){
        this->c = c;
    }
    
    void setNIter(unsigned int N_iter){
        this->N_iter = N_iter;
    }

    void setId(std::pair<int,int> id){
        this->id = id;
    }

    std::pair<int,int> getId(){
        return id;
    }

    bool getAns(){
        return ans;
    }

    void run() override {
        std::complex<float> zk(0.0f, 0.0f);
        for (int i = 1; i <= N_iter; i++) {
            zk = zk * zk + c;
            if (std::abs(zk) > 2.0) {
                ans = false;
                return;
            }
        }
        ans = true;
    }

    void sendTask(unsigned int to) override {
        MPI_Send(&N_iter, 1, MPI_UNSIGNED, to, MPI_TAG_MANDELBROT_TASK, MPI_COMM_WORLD);
        MPI_Send(&c, 1, MPI_C_COMPLEX, to, MPI_TAG_MANDELBROT_TASK, MPI_COMM_WORLD);
    }

    void getTask(unsigned int from) override {
        MPI_Recv(&N_iter, 1, MPI_UNSIGNED, from, MPI_TAG_MANDELBROT_TASK, MPI_COMM_WORLD, nullptr);
        MPI_Recv(&c, 1, MPI_C_COMPLEX, from, MPI_TAG_MANDELBROT_TASK, MPI_COMM_WORLD, nullptr);
    }

    void getResult(unsigned int from) override {
        MPI_Recv(&ans, 1, MPI_C_BOOL, from, MPI_TAG_MANDELBROT_TASK, MPI_COMM_WORLD, nullptr);
    }

    void sendResult(unsigned int to) override {
        MPI_Send(&ans, 1, MPI_C_BOOL, to, MPI_TAG_MANDELBROT_TASK, MPI_COMM_WORLD);
    }
};

class Mandelbrot {
public:
    Mandelbrot(std::pair<float, float> x_range, std::pair<float, float> y_range,
               unsigned int x_part, unsigned int y_part, unsigned int N_iter,
               unsigned int N_slaves)
        : x_range(x_range), y_range(y_range), x_part(x_part), y_part(y_part),
          N_iter(N_iter), N_slaves(N_slaves) {
        hasRun = false;
        fractal.resize(x_part, std::vector<bool>(y_part));
    }

    std::vector<std::vector<bool>> getFractal() {
        if (!hasRun)
            run();
        return fractal;
    }

private:
    std::pair<float, float> x_range, y_range;
    unsigned int x_part, y_part;
    unsigned int N_iter;
    bool hasRun;
    unsigned int N_slaves;

    std::vector<std::vector<bool>> fractal;

    void run() {
        MPI_Status status;
        std::pair<int,int> id;

        int pend_task = 0;
        int curr_slave = 1;

        hasRun = true;
        float x_delta = x_range.second - x_range.first;
        float y_delta = y_range.second - y_range.first;

        MandelbrotTask task[N_slaves + 1];
        for(int i=1; i<=N_slaves; i++)
            task[i].setNIter(N_iter);

        for (int i = 0; i < x_part; i++)
            for (int j = 0; j < y_part; j++) {
                float x = x_range.first + x_delta * ((float)i / (x_part - 1.0));
                float y = y_range.first + y_delta * ((float)j / (y_part - 1.0));
                

                if(curr_slave <= N_slaves){ // Assign a task to each slave
                    id = std::make_pair(i, j);
                    task[curr_slave].setC(std::complex<float>(x,y));
                    task[curr_slave].setId(id);
                    task[curr_slave].sendTask(curr_slave);
                    curr_slave++;
                    pend_task++;
                }else{ // Receive and send task
                    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
                    id = task[status.MPI_SOURCE].getId();
                    task[status.MPI_SOURCE].getResult(status.MPI_SOURCE);
                    fractal[id.first][id.second] = task[status.MPI_SOURCE].getAns();
                    pend_task--;

                    id = std::make_pair(i, j);
                    task[status.MPI_SOURCE].setC(std::complex<float>(x,y));
                    task[status.MPI_SOURCE].setId(id);
                    task[status.MPI_SOURCE].sendTask(status.MPI_SOURCE);
                    pend_task++;
                }
            }

        while(pend_task){
            MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            id = task[status.MPI_SOURCE].getId();
            task[status.MPI_SOURCE].getResult(status.MPI_SOURCE);
            fractal[id.first][id.second] = task[status.MPI_SOURCE].getAns();
            pend_task--;
        }

        for(int i=1; i<=N_slaves; i++)
            task[i].finalizeProcess(i);
    }
};

void exportToPGM(const std::vector<std::vector<bool>> &mat,
                 const std::string &filepath) {

    std::streambuf *cout_backup, *file_buf;
    std::ofstream filestr;

    filestr.open(filepath.c_str());
    cout_backup = std::cout.rdbuf();
    file_buf = filestr.rdbuf();
    std::cout.rdbuf(file_buf);

    int n = mat.size();
    int m = mat[0].size();

    std::cout << "P2" << std::endl;
    std::cout << n << ' ' << m << std::endl;
    std::cout << 1 << std::endl;

    for (int j = 0; j < m; j++)
        for (int i = 0; i < n; i++) {
            std::cout << mat[i][j];
            if (i + 1 < n)
                std::cout << " ";
            else
                std::cout << std::endl;
        }
    std::cout.rdbuf(cout_backup);

    filestr.close();
}

int mainMaster(int argc, char *argv[], int world_size) {
    auto start_time = std::chrono::high_resolution_clock::now();

    Mandelbrot F(std::make_pair(-1.5, 0.5), std::make_pair(-1.0, 1.0), 1024,
                 1024, 1'000'000, world_size - 1);
    auto mat = F.getFractal();

    auto end_time = std::chrono::high_resolution_clock::now();
    double total_time =
        std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time)
            .count();

    exportToPGM(mat, "fractal_par.pgm");

    freopen("output_par.txt", "w", stdout);
    std::cout << "Total time: " << total_time << " s" << std::endl;
    fclose(stdout);

    MPI_Finalize();

    return 0;
}

int mainSlave(int rank) {
    Task *task;
    while (true) {
        MPI_Status status;
        MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        if (status.MPI_TAG == MPI_TAG_END_PROCESS) {
            break;
        } else if (status.MPI_TAG == MPI_TAG_MANDELBROT_TASK)
            task = new MandelbrotTask;

        task->getTask(0);
        task->run();
        task->sendResult(0);
    }
    MPI_Finalize();
    return 0;
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_size < 2) {
        std::cerr << "Requires at least two processes." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
        exit(-1);
    }

    if (world_rank == 0) { // Master
        return mainMaster(argc, argv, world_size);
    } else {
        return mainSlave(world_rank);
    }
}