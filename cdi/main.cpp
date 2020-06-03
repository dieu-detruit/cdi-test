#include <fftw3.h>
#include <grid/bundle.hpp>
#include <grid/core.hpp>
#include <grid/linear.hpp>
#include <unit/double.hpp>

#include <omp.h>

#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

#include "constants.hpp"

template <class T>
bool is_in_closed(T x, T min, T max)
{
    return min <= x && x <= max;
}

template <class grid_array>
void array_print_file(grid_array& ar, std::string filename)
{
    std::ofstream file(filename);
    for (auto x : ar.template line<0>()) {
        for (auto y : ar.template line<1>()) {
            file << x / 1.0_m << ' '
                 << y / 1.0_m << ' '
                 << ar.at(x, y).arg() / 1.0_rad << ' '
                 << (ar.at(x, y) * ar.at(x, y).conj()).real() / amp_unit / amp_unit << std::endl;
        }
        file << std::endl;
    }
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
                 << (ar.at(x, y) * ar.at(x, y).conj()).real() / amp_unit / amp_unit << std::endl;
        }
        file << std::endl;
    }
}

int main()
{
    if (!fftw_init_threads()) {
        std::cerr << "error at initialization of threads" << std::endl;
        return 1;
    }
    fftw_plan_with_nthreads(4);

    std::cout << "making a problem..." << std::endl;

    Grid::DynamicRange<Length> detector_range{-0.5 * detector_length, 0.5 * detector_length, detector_pixel_num};
    Grid::DynamicRange<Length> diffraction_plane_range{-0.5 * diffraction_plane_length, 0.5 * diffraction_plane_length, detector_pixel_num};

    Grid::StaticRange<Length, detector_pixel_num> static_detector_range{
        -0.5 * detector_length, 0.5 * detector_length};
    Grid::StaticRange<Length, detector_pixel_num> static_diffraction_plane_range{
        -0.5 * diffraction_plane_length, 0.5 * diffraction_plane_length};

    std::random_device seed_gen{};
    std::mt19937 engine{seed_gen()};

    // 強度マップ(既知)
    Grid::GridVector<PhotonFluxDensity, Length, 2> detected_I{detector_range, detector_range};
    Grid::GridVector<WaveAmplitude, Length, 2> detected_I_sqrt{detector_range, detector_range};

    {  // 実際のマップたち(未知)

        // それなりに疎な出口関数
        Grid::GridVector<Complex<WaveAmplitude>, Length, 2> exit{diffraction_plane_range, diffraction_plane_range};
        {
            Grid::DynamicRange<Length> object_range{-0.5 * object_length, 0.5 * object_length, image_pixel_num};
            Grid::GridVector<Complex<WaveAmplitude>, Length, 2> image{object_range, object_range};
            std::ifstream file("../image/transobj.txt");

            for (auto& px : image) {
                int pixel;
                file >> pixel;  // pixel value is in [0, 255]
                px = 1.0 * pixel / 255.0 * amp_unit;
            }

            exit.fill(0.0 * amp_unit);

            for (auto [x, y] : image.lines()) {
                exit.at(x, y) += image.at(x, y);
            }
            print_file(exit, "exit.txt");
        }

        // FFT 出口関数 -> ディテクター面での波面
        Grid::GridVector<Complex<WaveAmplitude>, Length, 2> detected_known{detector_range, detector_range};
        {
            fftw_plan plan = fftw_plan_dft_2d(detector_pixel_num, detector_pixel_num,
                reinterpret_cast<fftw_complex*>(exit.data()), reinterpret_cast<fftw_complex*>(detected_known.data()),
                FFTW_FORWARD, FFTW_ESTIMATE);
            fftw_execute(plan);
            fftw_destroy_plan(plan);
            for (auto& d : detected_known) {
                d /= (double)(detector_pixel_num * detector_pixel_num);
            }
            print_file(detected_known, "detected_known.txt");
        }

        // 強度マップを取る {
        for (auto [f, I, rtI] : Grid::zip(detected_known, detected_I, detected_I_sqrt)) {
            I = f * f.conj() / (double)detector_pixel_num;
            rtI = I.sqrt();
        }
        // }

        std::ofstream problem_file("problem.txt");
        for (auto& x : detected_I.line(0)) {
            for (auto& y : detected_I.line(1)) {
                problem_file << x / 1.0_m << ' ' << y / 1.0_m << ' ' << detected_I.at(x, y) / amp_unit / amp_unit << std::endl;
            }
            problem_file << std::endl;
        }

        // もう一度戻して検証
        {
            fftw_plan plan = fftw_plan_dft_2d(detector_pixel_num, detector_pixel_num,
                reinterpret_cast<fftw_complex*>(detected_known.data()), reinterpret_cast<fftw_complex*>(exit.data()),
                FFTW_BACKWARD, FFTW_ESTIMATE);
            fftw_execute(plan);
            fftw_destroy_plan(plan);
            for (auto& e : exit) {
                e /= (double)(detector_pixel_num);
            }
            print_file(exit, "exit_test.txt");
        }
    }

    std::cout << "preparing to solve the problem..." << std::endl;

    // 以下で位相回復をする
    Grid::GridArray<Complex<WaveAmplitude>, Length, detector_pixel_num, detector_pixel_num> detected{static_detector_range, static_detector_range};
    Grid::GridArray<Complex<WaveAmplitude>, Length, detector_pixel_num, detector_pixel_num> exit{static_diffraction_plane_range, static_diffraction_plane_range};
    Grid::GridArray<Complex<WaveAmplitude>, Length, detector_pixel_num, detector_pixel_num> exit_prev{static_diffraction_plane_range, static_diffraction_plane_range};
    Grid::GridArray<bool, Length, detector_pixel_num, detector_pixel_num> support{static_diffraction_plane_range, static_diffraction_plane_range};

    // planの初期化は配列の初期化前にする(FFTW3公式より)
    fftw_plan plan_to_reciprocal
        = fftw_plan_dft_2d(detector_pixel_num, detector_pixel_num,
            reinterpret_cast<fftw_complex*>(exit.data()), reinterpret_cast<fftw_complex*>(detected.data()), FFTW_FORWARD, FFTW_MEASURE);
    fftw_plan plan_to_real
        = fftw_plan_dft_2d(detector_pixel_num, detector_pixel_num,
            reinterpret_cast<fftw_complex*>(detected.data()), reinterpret_cast<fftw_complex*>(exit.data()), FFTW_BACKWARD, FFTW_MEASURE);

    // サポートの設定
    for (auto [x, y] : support.lines()) {
        //support[x][y] = (x * x + y * y < 0.25 * object_diameter * object_diameter);
        support[x][y] = is_in_closed(x, -0.5 * object_length, 0.5 * object_length) and is_in_closed(y, -0.5 * object_length, 0.5 * object_length);
    }

    // 初期値を適当に撒く
    std::uniform_real_distribution<double> dist_init{0.0, 1.0};
    for (auto [in_support, e, e_prev] : Grid::zip(support, exit, exit_prev)) {
        e = e_prev = (in_support ? dist_init(engine) : 0.0) * amp_unit;
    }
    array_print_file(exit, "init.txt");

    exit.fill(0.0 * amp_unit);

    // 反復計算
    int i = 0;
    for (auto target : std::array<std::size_t, 6>{11, 21, 31, 41, 51, 101}) {
        for (; i < target; ++i) {
            if (i % 10 == 0) {
                std::cout << "rep: " << i << std::endl;
            }

            // copy to prev
            for (auto [e, e_prev] : Grid::zip(exit, exit_prev)) {
                e_prev = e;
            }

            // 出口関数 -> ディテクター面
            fftw_execute(plan_to_reciprocal);

            // 逆空間拘束
            for (auto [f, rtI] : Grid::zip(detected, detected_I_sqrt)) {
                Phase phase = f.arg();
                f = Complex<WaveAmplitude>::polar(rtI, phase);
            }

            // ディテクター面 -> 出口関数
            fftw_execute(plan_to_real);

            if (i == target - 1) {
                array_print_file(exit, "epoch_" + std::to_string(i) + "_exit.txt");
            }

            // 実空間拘束
            for (auto [in_support, e, e_prev] : Grid::zip(support, exit, exit_prev)) {
                e = in_support ? e : e_prev - beta * e;
            }
        }
    }

    fftw_destroy_plan(plan_to_real);
    fftw_destroy_plan(plan_to_reciprocal);

    // 最後にディテクター面から出口関数に戻す
    {
        fftw_plan plan = fftw_plan_dft_2d(detector_pixel_num, detector_pixel_num,
            reinterpret_cast<fftw_complex*>(detected.data()), reinterpret_cast<fftw_complex*>(exit.data()), FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);
    }

    array_print_file(exit, "answer.txt");

    fftw_cleanup_threads();

    return 0;
}
