#include <fftw3.h>

#include <array>
#include <filesystem>
#include <iostream>
#include <random>

#include <grid/core.hpp>
#include <grid/zip.hpp>
#include <unit/double.hpp>

int main()
{
    std::cout << std::endl;
    double itr_max = 100000000;
    std::cout << '[';
    for (int i = 0; i < 50; ++i) {
        std::cout << ' ';
    }
    std::cout << "]\r[";

    std::size_t hoge = itr_max / 50;
    int fuga;
    for (int i = 0; i < itr_max; ++i) {
        fuga = i * i;
        if (i % hoge == 0) {
            std::cout << '*';
        }
    }
    std::cout << std::endl;

    return 0;
}

#if 0
int main()
{
    using complex_t = std::complex<double>;

    std::cout << "sizeof(double) = " << sizeof(double) << std::endl;
    std::cout << "sizeof(complex_t) = " << sizeof(complex_t) << std::endl;

    complex_t* hoge = (complex_t*)fftw_malloc(sizeof(fftw_complex) * 20);
    for (std::size_t i = 0; i < 4; ++i) {
        for (std::size_t j = 0; j < 5; ++j) {
            std::cout << std::hex << &hoge[5 * i + j] << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "grid: " << std::endl;

    Grid::GridVector<complex_t, double, 2> fuga{{-1.0, 1.0, 4}, {-1.0, 1.0, 5}};
    for (auto& x : fuga.line(0)) {
        for (auto& y : fuga.line(1)) {
            std::cout << std::hex << &fuga.at(x, y) << ' ';
        }
        std::cout << std::endl;
    }

    return 0;
}
#endif

#if 0
int main(int argc, const char* argv[])
{
    using complex_t = std::complex<double>;

    constexpr std::size_t size_h = 32;
    constexpr std::size_t size_w = 32;
    const int fftsize = size_h * size_w;

    fftw_plan plan_fft, plan_ifft;
    int i;

    complex_t* in = (complex_t*)fftw_malloc(sizeof(complex_t) * fftsize);
    complex_t* out = (complex_t*)fftw_malloc(sizeof(complex_t) * fftsize);

    plan_fft = fftw_plan_dft_2d(32, 32, (fftw_complex*)in, (fftw_complex*)out, FFTW_FORWARD, FFTW_ESTIMATE);
    plan_ifft = fftw_plan_dft_2d(32, 32, (fftw_complex*)in, (fftw_complex*)out, FFTW_BACKWARD, FFTW_ESTIMATE);

    /* create sin wave of 440Hz */
    double sum_value_before = 0.0;
    for (auto [i, j] : Grid::zip(Grid::arange(0, 32, 1), Grid::arange(0, 32, 1))) {
        double x = 2.0 / 32 * (i - 16);
        double y = 2.0 / 32 * (j - 16);
        in[size_h * i + j] = std::exp(std::complex<double>(0.0, 440.0 / 44100.0 * 2.0 * M_PI * 32.0 * (x + y)));
        sum_value_before += std::norm(in[size_h * i + j]);
    }

    fftw_execute(plan_fft);
    fftw_execute(plan_ifft);

    /* get amplitude */
    double sum_value_after = 0.0;
    for (auto [i, j] : Grid::zip(Grid::arange(0, 32, 1), Grid::arange(0, 32, 1))) {
        sum_value_after += std::norm(in[size_h * i + j]);
    }

    std::cout << "before: " << sum_value_before << std::endl;
    std::cout << "after:  " << sum_value_after << std::endl;

    fftw_destroy_plan(plan_fft);
    fftw_destroy_plan(plan_ifft);
    fftw_free(in);
    fftw_free(out);
    return 0;
}
#endif

#if 0
int main(int argc, const char* argv[])
{
    using complex_t = std::complex<double>;

    constexpr std::size_t size_h = 32;
    constexpr std::size_t size_w = 32;
    const int fftsize = size_h * size_w;

    fftw_plan plan_fft, plan_ifft;
    int i;

    Grid::GridVector<complex_t, double, 2> in{{-1.0, 1.0, size_w}, {-1.0, 1.0, size_h}};
    Grid::GridVector<complex_t, double, 2> out{{-1.0, 1.0, size_w}, {-1.0, 1.0, size_h}};

    plan_fft = fftw_plan_dft_2d(32, 32, (fftw_complex*)in.data(), (fftw_complex*)out.data(), FFTW_FORWARD, FFTW_ESTIMATE);
    plan_ifft = fftw_plan_dft_2d(32, 32, (fftw_complex*)out.data(), (fftw_complex*)in.data(), FFTW_BACKWARD, FFTW_ESTIMATE);

    /* create sin wave of 440Hz */
    double sum_value_before = 0.0;
    for (auto [x, y] : in.lines()) {
        in[x][y] = 2.5 * std::exp(std::complex<double>(0.0, 440.0 / 44100.0 * 2.0 * M_PI * (x + y)));
        sum_value_before += std::norm(in[x][y]);
    }

    fftw_execute(plan_fft);
    fftw_execute(plan_ifft);

    double sum_value_after = 0.0;
    for (auto& v : in) {
        sum_value_after += std::norm(v);
    }

    fftw_execute(plan_fft);
    fftw_execute(plan_ifft);

    double sum_value_after_2 = 0.0;
    for (auto& v : in) {
        sum_value_after_2 += std::norm(v);
    }

    std::cout << "before: " << sum_value_before << std::endl;
    std::cout << "after / N^2:  " << sum_value_after / fftsize / fftsize << std::endl;
    std::cout << "after_2 / N^2:  " << sum_value_after / fftsize / fftsize << std::endl;


    fftw_destroy_plan(plan_fft);
    fftw_destroy_plan(plan_ifft);
    return 0;
}
#endif


#if 0
template <class T>
bool is_in_closed(T x, T min, T max)
{
    return min <= x && x <= max;
}

int main()
{
    using namespace Unit;

    Grid::GridVector<std::complex<double>, double, 2> in{{-1.0, 1.0, 256}, {-1.0, 1.0, 256}};
    Grid::GridVector<std::complex<double>, double, 2> out{{-1.0, 1.0, 256}, {-1.0, 1.0, 256}};

    fftw_plan plan_fft = fftw_plan_dft_2d(256, 256,
        reinterpret_cast<fftw_complex*>(in.data()), reinterpret_cast<fftw_complex*>(out.data()),
        FFTW_FORWARD, FFTW_ESTIMATE);

    fftw_plan plan_ifft = fftw_plan_dft_2d(256, 256,
        reinterpret_cast<fftw_complex*>(out.data()), reinterpret_cast<fftw_complex*>(in.data()),
        FFTW_BACKWARD, FFTW_ESTIMATE);

    std::random_device rd{};
    std::mt19937 gen{rd()};

    std::normal_distribution<double> dist_x{0.0, 0.3};
    std::normal_distribution<double> dist_y{0.0, 0.3};

    in.fill(0.0);
    out.fill(0.0);

    for (auto [x, y] : in.lines()) {
        in[x][y] = std::exp(std::complex<double>(0.0, 440.0 / 44100.0 * 2.0 * M_PI * (x + y)));
    }

    //constexpr std::size_t source_sample = 50000;
    //for (std::size_t i = 0; i < source_sample; ++i) {
    //double x = dist_x(gen), y = dist_y(gen);
    //if (is_in_closed(x, -1.0, 1.0) && is_in_closed(y, -1.0, 1.0)) {
    //in.at(x, y) += 1.0;
    //}
    //}
    double sum_value_before = 0.0;
    for (auto& v : in) {
        sum_value_before += std::norm(v);
    }
    std::cout << "before: " << sum_value_before << std::endl;

    std::cout << "calculating ..." << std::endl;

    fftw_execute(plan_fft);
    fftw_execute(plan_ifft);

    double sum_value_after = 0.0;
    for (auto& v : in) {
        sum_value_after += std::norm(v);
    }

    std::cout << "after:  " << sum_value_after << std::endl;

    fftw_destroy_plan(plan_fft);
    fftw_destroy_plan(plan_ifft);
    return 0;
}
#endif
