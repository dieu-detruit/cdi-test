#include <fftw3.h>
#include <grid/bundle.hpp>
#include <grid/core.hpp>
#include <grid/linear.hpp>
#include <stdexcept>
#include <unit/double.hpp>

#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>

#include "constants.hpp"

template <class T>
bool is_in_closed(T x, T min, T max)
{
    return min <= x && x <= max;
}

template <class complex_grid_vector>
void print_file(complex_grid_vector& ar, std::string filename)
{
    std::ofstream file("data_knifeedge/" + filename);
    for (auto x : ar.line(0)) {
        for (auto y : ar.line(1)) {
            file << x / 1.0_m << ' '
                 << y / 1.0_m << ' '
                 << std::arg(ar.at(x, y)) << ' ';
            if constexpr (std::is_same_v<decltype(ar.at(x, y)), std::complex<double>>) {
                file << std::norm(ar.at(x, y)) << std::endl;
            } else {
                file << std::norm(ar.at(x, y)).value << std::endl;
            }
        }
        file << std::endl;
    }
}

#define LOG_INFO(x) std::cout << #x ":" << (std::string(#x).length() < 15 ? "\t\t\t" : "\t\t") << x << std::endl;

struct DensityDistMeasurement {
    enum Direction {
        UP,
        DOWN,
        RIGHT,
        LEFT
    };

    static std::array<Direction, 4> dirs()
    {
        return {{UP, DOWN, RIGHT, LEFT}};
    }


    Direction dir;
    Length knife_pos;
    Grid::GridVector<WaveAmplitude, Length, 2> dist_abs;

    bool window(Length x, Length y)
    {
        switch (dir) {
        case UP:
            return y > knife_pos;
        case DOWN:
            return y < knife_pos;
        case RIGHT:
            return x > knife_pos;
        case LEFT:
            return x < knife_pos;
        }
        return false;
    };

    DensityDistMeasurement(Direction dir, Length knife_pos, Grid::DynamicRange<Length> range)
        : dir(dir), knife_pos(knife_pos), dist_abs{range, range} {}
};

int main()
{
    if (!fftw_init_threads()) {
        std::cerr << "error at initialization of threads" << std::endl;
        return 1;
    }
    fftw_plan_with_nthreads(4);

    std::cout << "making a problem..." << std::endl;

    std::filesystem::create_directory("data_knifeedge");

    Grid::DynamicRange<Length> detector_range{-0.5 * detector_length, 0.5 * detector_length, detector_pixel_num};
    Grid::DynamicRange<Length> probe_range{-0.5 * probe_length, 0.5 * probe_length, detector_pixel_num};

    std::random_device seed_gen{};
    std::mt19937 engine{seed_gen()};

    // 既知の情報
    Grid::GridVector<WaveAmplitude, Length, 2> dist_abs_open{detector_range, detector_range};
    std::vector<DensityDistMeasurement> detected_rtI_list;

    // 問題サンプルの生成
    {  // 実際のマップたち(未知)

        // 照明
        Grid::GridVector<Complex<WaveAmplitude>, Length, 2>
            src_image{{-0.5 * probe_length, 0.5 * probe_length, image_pixel_num}, {-0.5 * probe_length, 0.5 * probe_length, image_pixel_num}};
        {
            std::ifstream file("../image/transobj.txt");
            for (auto& px : src_image) {
                int pixel;
                file >> pixel;  // pixel value is in [0, 255]
                px = std::polar(pixel * amp_unit, 0.8_rad * M_PI * pixel / 255.0);
            }
        }
        print_file(src_image, "src_image.txt");

        Grid::GridVector<Complex<WaveAmplitude>, Length, 2> probe{probe_range, probe_range};
        probe.fill(0.0 * amp_unit);
        {
            for (auto [x, y] : probe.lines()) {
                probe.at(x, y) = src_image.at(x, y);
            }
        }

        print_file(probe, "probe.txt");

        // サンプル生成
        std::cout << "measurement..." << std::endl;
        {
            auto movement_range = Grid::linspace(-0.5 * probe_length + knife_step, 0.5 * probe_length - knife_step, knife_num - 1);
            Grid::GridVector<Complex<WaveAmplitude>, Length, 2> exit{probe_range, probe_range};
            Grid::GridVector<Complex<WaveAmplitude>, Length, 2> detected{detector_range, detector_range};

            auto exit_lines = exit.lines();
            fftw_plan plan = fftw_plan_dft_2d(detector_pixel_num, detector_pixel_num,
                (fftw_complex*)exit.data(), (fftw_complex*)detected.data(),
                FFTW_FORWARD, FFTW_ESTIMATE);
            int c = 0;

            // 画像全体
            {
                for (auto [x, y, e, p] : Grid::zip(exit_lines, exit, probe)) {
                    e = p;
                }

                // 出口 -> ディテクター面
                fftw_execute(plan);

                for (auto [f, rtI] : Grid::zip(detected, dist_abs_open)) {
                    rtI = std::abs(f) / (double)detector_pixel_num;
                }
                print_file(detected, "whole_measurement.txt");
            }

            // ナイフエッジあり
            for (auto dir : DensityDistMeasurement::dirs()) {
                for (auto& knife_pos : movement_range) {

                    const auto window = [&dir, &knife_pos](Length x, Length y) {
                        switch (dir) {
                        case DensityDistMeasurement::UP:
                            return y > knife_pos;
                        case DensityDistMeasurement::DOWN:
                            return y < knife_pos;
                        case DensityDistMeasurement::RIGHT:
                            return x > knife_pos;
                        case DensityDistMeasurement::LEFT:
                            return x < knife_pos;
                        }
                        return false;
                    };

                    // 出口波面 = 照射関数 * 窓
                    for (auto [x, y, e, p] : Grid::zip(exit_lines, exit, probe)) {
                        if (window(x, y)) {
                            e = p;
                        } else {
                            e = 0.0 * amp_unit;
                        }
                    }

                    if (c == 25) {
                        print_file(exit, "exit_example.txt");
                    }

                    // 出口 -> ディテクター面
                    fftw_execute(plan);

                    detected_rtI_list.emplace_back(dir, knife_pos, detector_range);

                    for (auto [f, rtI] : Grid::zip(detected, detected_rtI_list.back().dist_abs)) {
                        rtI = std::abs(f) / (double)detector_pixel_num;
                    }

                    if (c == 25) {
                        print_file(detected, "measurement.txt");
                    }
                    ++c;
                }
            }
            fftw_destroy_plan(plan);
            std::cout << "sample size: " << detected_rtI_list.size() << std::endl;
        }
    }


    std::cout << "Calculating ..." << std::endl;

    // タイコグラフィをやる
    Grid::GridVector<Complex<WaveAmplitude>, Length, 2> probe{probe_range, probe_range};

    // 反復計算
    Grid::GridVector<Complex<WaveAmplitude>, Length, 2> exit{probe_range, probe_range};
    Grid::GridVector<Complex<WaveAmplitude>, Length, 2> detected{detector_range, detector_range};

    auto probe_lines = probe.lines();

    fftw_plan plan_to_reciprocal = fftw_plan_dft_2d(probe_num, detector_pixel_num,
        (fftw_complex*)exit.data(), (fftw_complex*)detected.data(), FFTW_FORWARD, FFTW_MEASURE);
    fftw_plan plan_to_real = fftw_plan_dft_2d(detector_pixel_num, probe_num,
        (fftw_complex*)detected.data(), (fftw_complex*)exit.data(), FFTW_BACKWARD, FFTW_MEASURE);

    // プローブの初期化
    {
        Grid::GridVector<Complex<WaveAmplitude>, Length, 2> detected_open{detector_range, detector_range};
        std::uniform_real_distribution<double> dist_probe_init{0.0, 1.0};
        for (auto [d, rt_I] : Grid::zip(detected_open, dist_abs_open)) {
            d = rt_I / detector_pixel_num;
        }
        fftw_plan plan = fftw_plan_dft_2d(detector_pixel_num, probe_num,
            (fftw_complex*)detected_open.data(), (fftw_complex*)probe.data(), FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(plan);

        print_file(probe, "probe_init.txt");
        fftw_destroy_plan(plan);
    }

    bool first = true;

    for (int iteration = 0; iteration < 201; ++iteration) {

        std::cout << "Rep: " << iteration << std::endl;

        int i = 0;

        // ナイフなし
        {
            for (auto [p, e] : Grid::zip(probe, exit)) {
                e = p;
            }

            // 出口 -> ディテクター面 (exit -> detected)
            fftw_execute(plan_to_reciprocal);

            // 逆空間拘束
            constexpr double fft_scale_factor = detector_pixel_num * detector_pixel_num;
            for (auto [f, rtI] : Grid::zip(detected, dist_abs_open)) {
                f = std::polar(rtI / fft_scale_factor, f.arg());
            }

            // ディテクター面 -> 出口 (detected(constrained) -> exit)
            fftw_execute(plan_to_real);

            // 実空間での更新 (Probe)
            for (auto [p, e] : Grid::zip(probe, exit)) {
                p += 0.1 * (e - p);
            }
        }

        for (auto& mes : detected_rtI_list) {

            // ナイフを突き立てる
            for (auto [x, y, p, e] : Grid::zip(probe_lines, probe, exit)) {
                if (mes.window(x, y)) {
                    e = p;
                } else {
                    e = 0.0 * amp_unit;
                }
            }

            if (first) {
                print_file(exit, "itr_exit_example.txt");

                std::cout << "knife_pos: " << mes.knife_pos << std::endl;
                std::cout << "dir: " << mes.dir << std::endl;
                first = false;
            }

            // 出口 -> ディテクター面 (exit -> detected)
            fftw_execute(plan_to_reciprocal);

            // 逆空間拘束
            constexpr double fft_scale_factor = detector_pixel_num * detector_pixel_num;
            for (auto [f, rtI] : Grid::zip(detected, mes.dist_abs)) {
                f = std::polar(rtI / fft_scale_factor, f.arg());
            }

            // ディテクター面 -> 出口 (detected(constrained) -> exit)
            fftw_execute(plan_to_real);

            // 実空間での更新 (Probe)
            for (auto [x, y, p, e] : Grid::zip(probe_lines, probe, exit)) {
                if (mes.window(x, y)) {
                    p += 0.1 * (e - p);
                }
            }
        }

        if (iteration % 10 == 0) {
            print_file(probe, "epoch_" + std::to_string(iteration) + "_probe.txt");
        }
    }

    fftw_destroy_plan(plan_to_real);
    fftw_destroy_plan(plan_to_reciprocal);

    fftw_cleanup_threads();

    return 0;
}
