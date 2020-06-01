#include <fftw3.h>
#include <grid/core.hpp>
#include <grid/linear.hpp>
#include <grid/zip.hpp>

#include <cmath>
#include <complex>
#include <iostream>
#include <random>

bool is_in_closed(double x, double min, double max)
{
    return min <= x && x <= max;
}

using complex_t = std::complex<double>;

template <class T, std::size_t size1, std::size_t size2>
void output_array(const Grid::GridArray<T, double, size, size>& ar)
{
}

int main()
{
    constexpr std::size_t light_size = 50;
    constexpr std::size_t obj_size = 225;
    constexpr std::size_t step = 25;

    // オブジェクトを[-1,1]^2にマップします

    constexpr double size_ratio = (double)light_size / obj_size;

    constexpr double sigma_x = 0.4 * size_ratio;
    constexpr double sigma_y = 0.4 * size_ratio;

    // make light source {
    Grid::GridArray<complex_t, double, light_size, light_size> source{{-size_ratio, size_ratio}, {-size_ratio, size_ratio}};
    source.fill(0.0);

    std::random_device rd{};
    std::mt19937 gen{rd()};

    std::normal_distribution<double> dist_x{0.0, sigma_x};
    std::normal_distribution<double> dist_y{0.0, sigma_y};

    constexpr std::size_t source_sample = 50000;
    for (std::size_t i = 0; i < source_sample; ++i) {
        double x = dist_x(gen), y = dist_y(gen);
        if (is_in_closed(x, -size_ratio, size_ratio) && is_in_closed(y, -size_ratio, size_ratio)) {
            source.at(x, y) += 1.0;
        }
    }

    // }

    // 透過関数全体
    Grid::GridArray<complex_t, double, obj_size, obj_size> transparent{range, range};
    {
        std::ifstream file("../image/transobj.txt");

        for (auto& t : transparent) {
            int x;
            file >> x;  // x is value in [0, 255]
            t.real(1.0);
            t.imag(-x / 255.0 * 0.5);
        }
    }

    constexpr double dx = 2.0 / size;
    constexpr double dy = 2.0 / size;

    for (std::size_t offset_x : Grid::arange(0, obj_size - light_size, step)) {
        for (std::size_t offset_y = 0; offset_y + light_size <= obj_size; offset_y += step) {

            for (auto& x : Grid::arange(-1.0 + size_ratio * offset_x, -1.0 + size_ratio * (offset_x + 1), 2.0 / obj_size)) {
                for (auto& y : Grid::arange(-1.0 + size_ratio * offset_y, -1.0 + size_ratio * (offset_x + 1), 2.0 / obj_size)) {
                    std::cout << x << ' ' << y << ' ' << source.at(x, y).real() << ' ' << source.at(x, y).imag() << std::endl;
                }
                std::cout << std::endl;
            }

            Grid::GridArray::Range range{-1.0 + size_ratio * offset_x, -1.0 + size_ratio * (offset_x + 1)};

            // 透過関数
            Grid::GridArray<complex_t, double, size, size> transparent{range, range};

            for (auto& t : transparent) {
                int x;
                std::cin >> x;  // x is value in [0, 255]
                t.real(1.0);
                t.imag(-x / 255.0 * 0.5);
            }

            Grid::GridArray<complex_t, double, size, size> decayed{{-1.0, 1.0}, {-1.0, 1.0}};

            for (auto& x : Grid::arange(-1.0, 1.0, 2.0 / size)) {
                for (auto& y : Grid::arange(-1.0, 1.0, 2.0 / size)) {
                    decayed.at(x, y) = source.at(x, y) * transparent.at(x, y);
                }
            }

            Grid::GridArray<complex_t, double, size, size> detected{{-1.0, 1.0}, {-1.0, 1.0}};
            fftw_plan plan;

            plan = fftw_plan_dft_2d(size, size, reinterpret_cast<fftw_complex*>(decayed.data()), reinterpret_cast<fftw_complex*>(detected.data()), FFTW_FORWARD, FFTW_ESTIMATE);
            fftw_execute(plan);

            Grid::GridArray<double, double, size, size> detected_I{{-1.0, 1.0}, {-1.0, 1.0}};
            for (auto& x : Grid::arange(-1.0, 1.0, 2.0 / size)) {
                for (auto& y : Grid::arange(-1.0, 1.0, 2.0 / size)) {
                }
                std::cout << std::endl;
            }


            fftw_destroy_plan(plan);
        }
    }


#if 0
    for (auto& x : Grid::arange(-1.0, 1.0, 2.0 / size)) {
        for (auto& y : Grid::arange(-1.0, 1.0, 2.0 / size)) {
            std::cout << x << ' ' << y << ' ' << transparent.at(x, y).real() << ' ' << transparent.at(x, y).imag() << std::endl;
        }
        std::cout << std::endl;
    }
#endif


#if 0
    for (auto& x : Grid::arange(-1.0, 1.0, 2.0 / size)) {
        for (auto& y : Grid::arange(-1.0, 1.0, 2.0 / size)) {
            std::cout << x << ' ' << y << ' ' << decayed.at(x, y).real() << ' ' << decayed.at(x, y).imag() << std::endl;
        }
        std::cout << std::endl;
    }
#endif

#if 0
    for (auto& x : Grid::arange(-1.0, 1.0, 2.0 / size)) {
        for (auto& y : Grid::arange(-1.0, 1.0, 2.0 / size)) {

            auto& tmp = detected.at(x, y);

            std::cout << x << ' ' << y << ' ' << tmp.real() << ' ' << tmp.imag() << ' ' << (tmp * std::conj(tmp)).real() << ' ' << std::arg(tmp) << std::endl;
        }
        std::cout << std::endl;
    }
#endif
}
