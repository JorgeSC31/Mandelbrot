#include <complex>
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>
#include <chrono>

class Mandelbrot {
public:
    Mandelbrot(std::pair<float, float> x_range, std::pair<float, float> y_range,
               unsigned int x_part, unsigned int y_part, unsigned int N_iter)
        : x_range(x_range), y_range(y_range), x_part(x_part), y_part(y_part),
          N_iter(N_iter) {
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

    std::vector<std::vector<bool>> fractal;

    void run() {
        hasRun = true;
        float x_delta = x_range.second - x_range.first;
        float y_delta = y_range.second - y_range.first;
        float x, y;
        std::complex<float> c;

        for (int i = 0; i < x_part; i++)
            for (int j = 0; j < y_part; j++) {
                x = x_range.first + x_delta * ((float)i / (x_part - 1.0));
                y = y_range.first + y_delta * ((float)j / (y_part - 1.0));

                c.real(x);
                c.imag(y);
                fractal[i][j] = inFractal(c);
            }
    }

    bool inFractal(std::complex<float> c) {
        std::complex<float> zk(0.0f, 0.0f);
        for (int i = 1; i <= N_iter; i++) {
            zk = zk * zk + c;
            if (std::abs(zk) > 2.0)
                return false;
        }
        return true;
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

int main(int argc, char *argv[]) {
    auto start_time = std::chrono::high_resolution_clock::now();

    Mandelbrot F(std::make_pair(-1.5, 0.5), std::make_pair(-1.0, 1.0), 1024, 1024,
                 1'000'000);
    auto mat = F.getFractal();
    
    auto end_time = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();

    exportToPGM(mat, "fractal_sec.pgm");
    freopen("output_sec.txt", "w", stdout);
    std::cout << "Total time: " << total_time << " s" << std::endl;
    fclose(stdout);
    return 0;
}