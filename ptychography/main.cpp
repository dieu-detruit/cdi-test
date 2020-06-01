#include <fftw3.h>
#include <grid/bundle.hpp>
#include <grid/core.hpp>
#include <grid/linear.hpp>
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

template <class grid_vector>
void print_file(grid_vector& ar, std::string filename)
{
    std::ofstream file(filename);
    for (auto x : ar.line(0)) {
        for (auto y : ar.line(1)) {
            file << x / 1.0_m << ' '
                 << y / 1.0_m << ' '
                 << ar.at(x, y).arg() / 1.0_rad << ' '
                 << (ar.at(x, y) * ar.at(x, y).conj()).real() / dens_unit << std::endl;
        }
        file << std::endl;
    }
}

#define LOG_INFO(x) std::cout << #x ":" << (std::string(#x).length() < 15 ? "\t\t\t" : "\t\t") << x << std::endl;

struct DensityDistMeasurement {
    Length offset_x;
    Length offset_y;
    Grid::GridVector<PhotonFluxDensity, Length, 2> dist;

    DensityDistMeasurement() : offset_x(0.0), offset_y(0.0),
                               dist{{-0.5 * detector_length, 0.5 * detector_length, detector_pixel_num}, {-0.5 * detector_length, 0.5 * detector_length, detector_pixel_num}} {}
    DensityDistMeasurement(Length offset_x, Length offset_y, Grid::DynamicRange<Length> range)
        : offset_x(offset_x), offset_y(offset_y), dist{range, range} {}
};

int main()
{
    if (!fftw_init_threads()) {
        std::cerr << "error at initialization of threads" << std::endl;
        return 1;
    }
    fftw_plan_with_nthreads(4);

    LOG_INFO(image_pixel_num);
    LOG_INFO(object_pixel_num);
    LOG_INFO(detector_pixel_num);
    LOG_INFO(light_pixel_num);

    LOG_INFO(detector_pixel_size);
    LOG_INFO(detector_length);
    LOG_INFO(object_length);
    LOG_INFO(object_pixel_size);
    LOG_INFO(light_length);
    LOG_INFO(light_stddev);

    LOG_INFO(overlap_ratio);
    LOG_INFO(scanning_step);
    LOG_INFO(image_num);
    LOG_INFO(scanning_half_length);

    std::cout << "making a problem..." << std::endl;

    std::filesystem::create_directory("data_ptychography");

    Grid::DynamicRange<Length> detector_range{-0.5 * detector_length, 0.5 * detector_length, detector_pixel_num};
    Grid::DynamicRange<Length> object_range{-0.5 * object_length, 0.5 * object_length, object_pixel_num};
    Grid::DynamicRange<Length> light_range{-0.5 * light_length, 0.5 * light_length, detector_pixel_num};

    std::random_device seed_gen{};
    std::mt19937 engine{seed_gen()};

    // 強度マップ(既知)
    std::vector<DensityDistMeasurement> detected_I_list;
    {  // 実際のマップたち(未知)

        // ガウシアン照明関数
        Grid::GridVector<Complex<WaveAmplitude>, Length, 2> light{light_range, light_range};
        std::normal_distribution<> dist_x(0.0, light_stddev / 1.0_m);
        std::normal_distribution<> dist_y(0.0, light_stddev / 1.0_m);

        constexpr std::size_t source_sample = 50000;
        for (std::size_t i = 0; i < source_sample; ++i) {
            Length x = dist_x(engine) * 1.0_m;
            Length y = dist_y(engine) * 1.0_m;
            if (is_in_closed(x, -0.5 * light_length, 0.5 * light_length) && is_in_closed(y, -0.5 * light_length, 0.5 * light_length)) {
                light.at(x, y) += 1.0 * amp_unit;
            }
        }

        // 疑似オブジェクト
        Grid::GridVector<std::complex<double>, Length, 2> src_image{{-0.5 * object_length, 0.5 * object_length, image_pixel_num}, {-0.5 * object_length, 0.5 * object_length, image_pixel_num}};
        {
            std::ifstream file("../image/transobj.txt");
            for (auto& px : src_image) {
                int pixel;
                file >> pixel;  // pixel value is in [0, 255]
                px = std::polar(0.95, 0.1 * M_PI * pixel / 255.0);
            }
        }

        Grid::GridVector<std::complex<double>, Length, 2> object{object_range, object_range};
        object.fill(0.0);
        for (auto [x, y] : src_image.lines()) {
            object.at(x, y) += src_image.at(x, y);
        }
        {
            std::ofstream file("data_ptychography/object.txt");
            for (auto x : object.line(0)) {
                for (auto y : object.line(1)) {
                    file << x / 1.0_m << ' '
                         << y / 1.0_m << ' '
                         << std::arg(object.at(x, y)) << ' '
                         << std::norm(object.at(x, y)) << std::endl;
                }
                file << std::endl;
            }
        }

        // サンプル生成
        {
            auto offset_range = Grid::arange(-scanning_half_length, scanning_half_length, scanning_step);
            Grid::GridVector<Complex<WaveAmplitude>, Length, 2> exit{light_range, light_range};
            Grid::GridVector<Complex<WaveAmplitude>, Length, 2> detected{detector_range, detector_range};
            fftw_plan plan = fftw_plan_dft_2d(detector_pixel_num, detector_pixel_num,
                reinterpret_cast<fftw_complex*>(exit.data()), reinterpret_cast<fftw_complex*>(detected.data()),
                FFTW_FORWARD, FFTW_MEASURE);

            for (auto [offset_x, offset_y] : Grid::prod(offset_range, offset_range)) {

                // 出口波面 = 照射関数 * 透過関数
                for (auto [x, y] : exit.lines()) {
                    exit.at(x, y) = light.at(x, y) * object.at(x + offset_x, y + offset_y);
                }

                // 出口 -> ディテクター面
                fftw_execute(plan);

                detected_I_list.emplace_back(offset_x, offset_y, Grid::DynamicRange<Length>{-scanning_half_length, scanning_half_length, image_num});

                for (auto [f, I] : Grid::zip(detected, detected_I_list.back().dist)) {
                    I = (f * f.conj()).real();
                }
            }
            fftw_destroy_plan(plan);
            std::cout << "sample size: " << detected_I_list.size() << std::endl;
        }
    }

    std::cout << "Calculating ..." << std::endl;

    // タイコグラフィをやる
    Grid::GridVector<Complex<WaveAmplitude>, Length, 2> light{light_range, light_range};
    Grid::GridVector<std::complex<double>, Length, 2> object{object_range, object_range};

    // 照明関数の初期値散布
    std::uniform_real_distribution<double> dist_light_init{0.0, 1.0};
    for (auto& l : light) {
        l = dist_light_init(engine) * amp_unit;
    }

    // オブジェクトの初期値散布
    std::uniform_real_distribution<double> dist_obj_init{0.0, 1.0};
    for (auto& o : object) {
        o = dist_obj_init(engine) * amp_unit;
    }

    // 反復計算
    Grid::GridVector<std::complex<double>, Length, 2> object_prev{object_range, object_range};

    Grid::GridVector<Complex<WaveAmplitude>, Length, 2> exit{light_range, light_range};
    Grid::GridVector<Complex<WaveAmplitude>, Length, 2> detected{detector_range, detector_range};

    Grid::GridVector<Complex<WaveAmplitude>, Length, 2> exit_prev{light_range, light_range};
    Grid::GridVector<Complex<WaveAmplitude>, Length, 2> light_prev{light_range, light_range};

    fftw_plan plan_to_real = fftw_plan_dft_2d(detector_pixel_num, detector_pixel_num,
        reinterpret_cast<fftw_complex*>(detected.data()), reinterpret_cast<fftw_complex*>(exit.data()),
        FFTW_FORWARD, FFTW_MEASURE);
    fftw_plan plan_to_reciprocal = fftw_plan_dft_2d(detector_pixel_num, detector_pixel_num,
        reinterpret_cast<fftw_complex*>(exit.data()), reinterpret_cast<fftw_complex*>(detected.data()),
        FFTW_BACKWARD, FFTW_MEASURE);

    for (int iteration = 0; iteration < 100; ++iteration) {

        double itr_size = detected_I_list.size();
        std::cout << "Rep: " << iteration << std::endl;

        for (auto& mes : detected_I_list) {
            // 出口波面 = 照射関数 * 透過関数
            for (auto [x, y] : exit.lines()) {
                exit.at(x, y) = light.at(x, y) * object.at(x + mes.offset_x, y + mes.offset_y);
                object_prev.at(x, y) = object.at(x + mes.offset_x, y + mes.offset_y);
            }
            for (auto [l, l_prev] : Grid::zip(light, light_prev)) {
                l_prev = l;
            }
            for (auto [e, e_prev] : Grid::zip(exit, exit_prev)) {
                e_prev = e;
            }

            // 出口 -> ディテクター面
            fftw_execute(plan_to_reciprocal);

            // 逆空間拘束
            for (auto [f, I] : Grid::zip(detected, mes.dist)) {
                f = polar(I.sqrt(), f.arg());
            }

            // ディテクター面 -> 出口
            fftw_execute(plan_to_real);

            // 絶対値の最大値を取る
            PhotonFluxDensity l_prev_normsq_max = 0.0 * dens_unit;
            for (auto& l_prev : light_prev) {
                auto norm_sq = (l_prev * l_prev.conj()).real();
                if (l_prev_normsq_max < norm_sq) {
                    l_prev_normsq_max = norm_sq;
                }
            }
            double o_prev_normsq_max = 0.0;
            for (auto& o_prev : object_prev) {
                auto norm_sq = std::real(o_prev * std::conj(o_prev));
                if (o_prev_normsq_max < norm_sq) {
                    o_prev_normsq_max = norm_sq;
                }
            }

            // 実空間での更新
            for (auto [x, y] : exit.lines()) {
                auto tmp = exit.at(x, y) - exit_prev.at(x, y);
                auto& l_prev = light_prev.at(x, y);
                auto& o_prev = object_prev.at(x, y);

                object.at(x + mes.offset_x, y + mes.offset_y)
                    = o_prev + alpha * std::complex<double>(l_prev.conj() / l_prev_normsq_max * tmp);
                light.at(x, y)
                    = l_prev + beta * std::conj(o_prev) / o_prev_normsq_max * tmp;
            }
        }

        std::ofstream file("data_ptychography/epoch_" + std::to_string(iteration) + "_object.txt");
        for (auto x : object.line(0)) {
            for (auto y : object.line(1)) {
                file << x / 1.0_m << ' '
                     << y / 1.0_m << ' '
                     << std::arg(object.at(x, y)) << ' '
                     << std::norm(object.at(x, y)) << std::endl;
            }
            file << std::endl;
        }
        print_file(light, "data_ptychography/epoch_" + std::to_string(iteration) + "_light.txt");
    }

    fftw_destroy_plan(plan_to_real);
    fftw_destroy_plan(plan_to_reciprocal);

    fftw_cleanup_threads();

    return 0;
}
