#include <fftw3.h>
#include <grid/core.hpp>
#include <grid/linear.hpp>

#include <cmath>
#include <complex>
#include <iostream>
#include <random>

using complex_t = std::complex<double>;

int main()
{
    constexpr std::size_t size = 225;

    constexpr double dx = 2.0 / size;
    constexpr double dy = 2.0 / size;

    constexpr double sigma_x = 0.4;
    constexpr double sigma_y = 0.4;

    Grid::GridArray<complex_t, double, size, size> source{{-1.0, 1.0}, {-1.0, 1.0}};
    source.fill(0.0);

    std::random_device rd{};
    std::mt19937 gen{rd()};

    std::normal_distribution<double> dist_x{0.0, sigma_x};
    std::normal_distribution<double> dist_y{0.0, sigma_y};

    constexpr std::size_t source_sample = 50000;
    for (std::size_t i = 0; i < source_sample; ++i) {
        double x = dist_x(gen), y = dist_y(gen);
        if (is_in_closed(x, -1.0, 1.0) && is_in_closed(y, -1.0, 1.0)) {
            source.at(x, y) += 1.0;
        }
    }
#if 0 
        for (auto& x : Grid::arange(-1.0, 1.0, 2.0 / size)) {
            for (auto& y : Grid::arange(-1.0, 1.0, 2.0 / size)) {
                std::cout << x << ' ' << y << ' ' << source.at(x, y).real() << ' ' << source.at(x, y).imag() << std::endl;
            }
            std::cout << std::endl;
        }
#endif

    Grid::GridArray<complex_t, double, size, size> transparent{{-1.0, 1.0}, {-1.0, 1.0}};

    for (auto& t : transparent) {
        int x;
        std::cin >> x;  // x is value in [0, 255]
        t.real(1.0);
        t.imag(-x / 255.0 * 0.5);
    }

#if 0
    for (auto& x : Grid::arange(-1.0, 1.0, 2.0 / size)) {
        for (auto& y : Grid::arange(-1.0, 1.0, 2.0 / size)) {
            std::cout << x << ' ' << y << ' ' << transparent.at(x, y).real() << ' ' << transparent.at(x, y).imag() << std::endl;
        }
        std::cout << std::endl;
    }
#endif

    Grid::GridArray<complex_t, double, size, size> decayed{{-1.0, 1.0}, {-1.0, 1.0}};

    for (auto& x : Grid::arange(-1.0, 1.0, 2.0 / size)) {
        for (auto& y : Grid::arange(-1.0, 1.0, 2.0 / size)) {
            decayed.at(x, y) = source.at(x, y) * transparent.at(x, y);
        }
    }

#if 0
    for (auto& x : Grid::arange(-1.0, 1.0, 2.0 / size)) {
        for (auto& y : Grid::arange(-1.0, 1.0, 2.0 / size)) {
            std::cout << x << ' ' << y << ' ' << decayed.at(x, y).real() << ' ' << decayed.at(x, y).imag() << std::endl;
        }
        std::cout << std::endl;
    }
#endif

    Grid::GridArray<complex_t, double, size, size> detected{{-1.0, 1.0}, {-1.0, 1.0}};
    fftw_plan plan;

    plan = fftw_plan_dft_2d(size, size, reinterpret_cast<fftw_complex*>(decayed.data()), reinterpret_cast<fftw_complex*>(detected.data()), FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

#if 0
    for (auto& x : Grid::arange(-1.0, 1.0, 2.0 / size)) {
        for (auto& y : Grid::arange(-1.0, 1.0, 2.0 / size)) {

            auto& tmp = detected.at(x, y);

            std::cout << x << ' ' << y << ' ' << tmp.real() << ' ' << tmp.imag() << ' ' << (tmp * std::conj(tmp)).real() << ' ' << std::arg(tmp) << std::endl;
        }
        std::cout << std::endl;
    }
#endif

    Grid::GridArray<double, double, size, size> detected_I{{-1.0, 1.0}, {-1.0, 1.0}};
    for (auto& x : Grid::arange(-1.0, 1.0, 2.0 / size)) {
        for (auto& y : Grid::arange(-1.0, 1.0, 2.0 / size)) {
        }
        std::cout << std::endl;
    }


    fftw_destroy_plan(plan);
}
