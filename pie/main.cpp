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
    std::ofstream file("data_pie/" + filename);
    for (auto x : ar.line(0)) {
        for (auto y : ar.line(1)) {
            file << x / 1.0_m << ' '
                 << y / 1.0_m << ' '
                 << std::arg(ar.at(x, y)) << ' '
                 << std::norm(ar.at(x, y)) << std::endl;
        }
        file << std::endl;
    }
}

#define LOG_INFO(x) std::cout << #x ":" << (std::string(#x).length() < 15 ? "\t\t\t" : "\t\t") << x << std::endl;

struct DensityDistMeasurement {
    Length offset_x;
    Length offset_y;
    Grid::GridVector<WaveAmplitude, Length, 2> dist_abs;

    DensityDistMeasurement() : offset_x(0.0), offset_y(0.0),
                               dist_abs{{-0.5 * detector_length, 0.5 * detector_length, detector_pixel_num}, {-0.5 * detector_length, 0.5 * detector_length, detector_pixel_num}} {}
    DensityDistMeasurement(Length offset_x, Length offset_y, Grid::DynamicRange<Length> range)
        : offset_x(offset_x), offset_y(offset_y), dist_abs{range, range} {}
};

int main()
{
    if (!fftw_init_threads()) {
        std::cerr << "error at initialization of threads" << std::endl;
        return 1;
    }
    fftw_plan_with_nthreads(4);

    LOG_INFO(image_pixel_num);
    LOG_INFO(scanning_object_pixel_num);
    LOG_INFO(detector_pixel_num);
    LOG_INFO(reconstruct_box_pixel_num);

    LOG_INFO(detector_pixel_size);
    LOG_INFO(detector_length);
    LOG_INFO(object_length);
    LOG_INFO(reconstruct_box_length);
    LOG_INFO(light_stddev);

    LOG_INFO(overlap_ratio);
    LOG_INFO(scanning_step);
    LOG_INFO(measure_num);
    LOG_INFO(scanning_half_length);

    std::cout << "making a problem..." << std::endl;

    std::filesystem::create_directory("data_pie");

    Grid::DynamicRange<Length> detector_range{-0.5 * detector_length, 0.5 * detector_length, detector_pixel_num};
    Grid::DynamicRange<Length> object_range{-0.5 * object_length, 0.5 * object_length, scanning_object_pixel_num};
    Grid::DynamicRange<Length> reconstruct_box_range{-0.5 * reconstruct_box_length, 0.5 * reconstruct_box_length, reconstruct_box_pixel_num};

    std::random_device seed_gen{};
    std::mt19937 engine{seed_gen()};

    // 既知の情報
    std::vector<DensityDistMeasurement> detected_rtI_list;
    Grid::GridVector<Complex<WaveAmplitude>, Length, 2> probe{reconstruct_box_range, reconstruct_box_range};

    // 問題サンプルの生成
    {  // 実際のマップたち(未知)

        // ガウシアン照明関数
        std::cout << "probe..." << std::endl;
        std::normal_distribution<> dist_x(0.0, light_stddev / 1.0_m);
        std::normal_distribution<> dist_y(0.0, light_stddev / 1.0_m);

        constexpr std::size_t source_sample = 10 * reconstruct_box_pixel_num * reconstruct_box_pixel_num;
        for (std::size_t i = 0; i < source_sample; ++i) {
            Length x = dist_x(engine) * 1.0_m;
            Length y = dist_y(engine) * 1.0_m;
            if (is_in_closed(x, -0.5 * reconstruct_box_length, 0.5 * reconstruct_box_length) && is_in_closed(y, -0.5 * reconstruct_box_length, 0.5 * reconstruct_box_length)) {
                probe.at(x, y) += 1.0 * amp_unit;
            }
        }
        print_file(probe, "probe.txt");

        // 疑似オブジェクト
        std::cout << "object..." << std::endl;
        Grid::GridVector<std::complex<double>, Length, 2>
            src_image{{-0.5 * object_length, 0.5 * object_length, image_pixel_num}, {-0.5 * object_length, 0.5 * object_length, image_pixel_num}};
        {
            std::ifstream file("../image/transobj.txt");
            for (auto& px : src_image) {
                int pixel;
                file >> pixel;  // pixel value is in [0, 255]
                px = std::polar(0.95, 0.8 * M_PI * pixel / 255.0);
            }
        }
        print_file(src_image, "src_image.txt");

        Grid::GridVector<std::complex<double>, Length, 2> object{object_range, object_range};
        object.fill(0.0);
        {
            for (auto [x, y] : object.lines()) {
                object.at(x, y) = src_image.at(x, y);
            }
        }

        print_file(object, "object.txt");

        // サンプル生成
        std::cout << "measurement..." << std::endl;
        {
            auto offset_range = Grid::linspace(-scanning_half_length, scanning_half_length, measure_num);
            Grid::GridVector<Complex<WaveAmplitude>, Length, 2> exit{reconstruct_box_range, reconstruct_box_range};
            Grid::GridVector<Complex<WaveAmplitude>, Length, 2> detected{detector_range, detector_range};

            auto exit_lines = exit.lines();
            fftw_plan plan = fftw_plan_dft_2d(detector_pixel_num, detector_pixel_num,
                (fftw_complex*)exit.data(), (fftw_complex*)detected.data(),
                FFTW_FORWARD, FFTW_MEASURE);
            bool first = true;

            Grid::GridVector<Complex<WaveAmplitude>, Length, 2> whiteboard{object_range, object_range};

            for (auto [offset_x, offset_y] : Grid::prod(offset_range, offset_range)) {

                {
                    auto range = Grid::arange(-5.0_nm, 5.0_nm, 1.0_nm);
                    for (auto [x, y] : Grid::zip(Grid::prod(range, range))) {
                        whiteboard.at(x + offset_x, y + offset_y) += 1000.0 * amp_unit;
                    }
                    for (auto [x, y] : Grid::prod(probe.line(0), probe.line(1))) {
                        whiteboard.at(x + offset_x, y + offset_y) += 100.0 * amp_unit;
                    }
                }

                // 出口波面 = 照射関数 * 透過関数
                for (auto [x, y, e, p] : Grid::zip(exit_lines, exit, probe)) {
                    e = p * object.at(x + offset_x, y + offset_y);
                }

                if (first) {
                    print_file(exit, "measurement_exit_example.txt");
                }

                // 出口 -> ディテクター面
                fftw_execute(plan);

                if (first) {
                    print_file(detected, "measurement.txt");
                }

                detected_rtI_list.emplace_back(offset_x, offset_y, detector_range);

                for (auto [f, rtI] : Grid::zip(detected, detected_rtI_list.back().dist_abs)) {
                    rtI = std::abs(f) / (double)reconstruct_box_pixel_num;
                }

                if (first) {
                    std::ofstream file("data_pie/measurement_rtI.txt");
                    for (auto x : detected_rtI_list.back().dist_abs.line(0)) {
                        for (auto y : detected_rtI_list.back().dist_abs.line(1)) {
                            file << x / 1.0_m << ' '
                                 << y / 1.0_m << ' '
                                 << detected_rtI_list.back().dist_abs.at(x, y) << std::endl;
                        }
                        file << std::endl;
                    }
                    first = false;
                }
            }
            fftw_destroy_plan(plan);
            print_file(whiteboard, "whiteboard.txt");
            std::cout << "sample size: " << detected_rtI_list.size() << std::endl;
        }
    }

    std::cout << "Calculating ..." << std::endl;

    // タイコグラフィをやる
    Grid::GridVector<std::complex<double>, Length, 2> object{object_range, object_range};

    // 反復計算
    Grid::GridVector<std::complex<double>, Length, 2> reconstruct_box{reconstruct_box_range, reconstruct_box_range};

    Grid::GridVector<Complex<WaveAmplitude>, Length, 2> exit{reconstruct_box_range, reconstruct_box_range};
    Grid::GridVector<Complex<WaveAmplitude>, Length, 2> exit_new{reconstruct_box_range, reconstruct_box_range};
    Grid::GridVector<Complex<WaveAmplitude>, Length, 2> detected{detector_range, detector_range};

    Grid::GridVector<ObjectUpdateWeight, Length, 2> object_update_weight{reconstruct_box_range, reconstruct_box_range};
    Grid::GridVector<double, Length, 2> reconstruct_box_norm{reconstruct_box_range, reconstruct_box_range};

    auto reconstruct_box_lines = reconstruct_box.lines();

    fftw_plan plan_to_reciprocal = fftw_plan_dft_2d(reconstruct_box_pixel_num, detector_pixel_num,
        (fftw_complex*)exit.data(), (fftw_complex*)detected.data(), FFTW_FORWARD, FFTW_MEASURE);
    fftw_plan plan_to_real = fftw_plan_dft_2d(detector_pixel_num, reconstruct_box_pixel_num,
        (fftw_complex*)detected.data(), (fftw_complex*)exit_new.data(), FFTW_BACKWARD, FFTW_MEASURE);

    // オブジェクトの初期値散布
    std::uniform_real_distribution<double> dist_obj_init{0.0, 1.0};
    for (auto& o : object) {
        o = std::polar(0.95, dist_obj_init(engine));
    }
    print_file(object, "object_init.txt");

    // 実空間での更新係数
    {
        Grid::GridVector<PhotonFluxDensity, Length, 2> probe_norm{reconstruct_box_range, reconstruct_box_range};
        PhotonFluxDensity p_norm_max = 0.0 * dens_unit;
        for (auto [p, p_norm] : Grid::zip(probe, probe_norm)) {
            p_norm = std::norm(p);
            if (p_norm > p_norm_max) {
                p_norm_max = p_norm;
            }
        }
        for (auto [p, p_norm, w] : Grid::zip(probe, probe_norm, object_update_weight)) {
            w = (p_norm / p_norm_max).sqrt() * std::conj(p) / (p_norm + epsilon);
        }
    }
    print_file(object_update_weight, "object_update_weight.txt");

    bool first = true;
    for (int iteration = 0; iteration < 100; ++iteration) {

        std::cout << "Rep: " << iteration << std::endl;

        int i = 0;
        for (auto& mes : detected_rtI_list) {

            // オブジェクトの切り出し
            for (auto [x, y, r] : Grid::zip(reconstruct_box_lines, reconstruct_box)) {
                r = object.at(x + mes.offset_x, y + mes.offset_y);
            }


            // 出口波面 = 照射関数 * 透過関数
            for (auto [p, r, e] : Grid::zip(probe, reconstruct_box, exit)) {
                e = p * r;
            }

            if (first) {
                print_file(reconstruct_box, "reconstruct_box_example.txt");
                print_file(exit, "exit_example.txt");
            }

            // 出口 -> ディテクター面 (exit -> detected)
            fftw_execute(plan_to_reciprocal);

            if (first) {
                print_file(detected, "reciprocal_example.txt");
            }

            // 逆空間拘束
            constexpr double fft_scale_factor = reconstruct_box_pixel_num * reconstruct_box_pixel_num;
            for (auto [f, rtI] : Grid::zip(detected, mes.dist_abs)) {
                f = std::polar(rtI / fft_scale_factor, f.arg());
            }

            if (first) {
                print_file(detected, "reciprocal_constrained_example.txt");
            }

            // ディテクター面 -> 出口 (detected(constrained) -> exit_new)
            fftw_execute(plan_to_real);

            if (first) {
                print_file(exit_new, "exit_new_example.txt");
            }

            if (first) {
                Grid::GridVector<Complex<WaveAmplitude>, Length, 2> exit_diff{reconstruct_box_range, reconstruct_box_range};
                for (auto [e, e_new, e_diff] : Grid::zip(exit, exit_new, exit_diff)) {
                    e_diff = e_new - e;
                }
                print_file(exit_diff, "exit_diff_example.txt");
                first = false;
            }

            // 実空間での更新(Object)
            for (auto [r, w, e, e_new] : Grid::zip(reconstruct_box, object_update_weight, exit, exit_new)) {
                r += std::complex<double>(w * (e_new - e));
            }

            // 切り出したオブジェクトを戻す
            for (auto [x, y, r] : Grid::zip(reconstruct_box_lines, reconstruct_box)) {
                object.at(x + mes.offset_x, y + mes.offset_y) = r;
            }
        }

        // print_file(probe, "epoch_" + std::to_string(iteration) + "_probe.txt");
        print_file(object, "epoch_" + std::to_string(iteration) + "_object.txt");
    }

    fftw_destroy_plan(plan_to_real);
    fftw_destroy_plan(plan_to_reciprocal);

    fftw_cleanup_threads();

    return 0;
}
