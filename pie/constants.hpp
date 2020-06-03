#pragma once

#include <unit/double.hpp>

inline namespace Constant
{


using namespace Unit;

using WaveNumber = decltype(1.0 / Length{});
using PhotonFluxDensity = Amount;
using WaveAmplitude = decltype(PhotonFluxDensity{}.sqrt());
using ObjectUpdateWeight = decltype(Complex<DimensionLessType>{} / WaveAmplitude{});

inline constexpr WaveAmplitude amp_unit{1.0};
inline constexpr PhotonFluxDensity dens_unit{1.0};

inline constexpr std::size_t image_pixel_num = 847;
inline constexpr std::size_t detector_pixel_num = 201;
inline constexpr std::size_t reconstruct_box_pixel_num = detector_pixel_num;

inline constexpr Length detector_pixel_size = 50.0_um;
inline constexpr Length detector_length = double(detector_pixel_num) * detector_pixel_size;
inline constexpr Length object_length = 1.9_um;

inline constexpr Length light_stddev = 200.0_nm;
inline constexpr Length reconstruct_box_length = 500.0_nm;
inline constexpr Length reconstruct_box_pixel_size = reconstruct_box_length / (double)reconstruct_box_pixel_num;

inline constexpr Length scanning_step = 100.0_nm;
inline constexpr Length scanning_half_length = 0.5 * (object_length - reconstruct_box_length);
inline constexpr double overlap_ratio = 1.0 - (double)(scanning_step / reconstruct_box_length);

//inline constexpr std::size_t measure_num = 1UL + std::size_t(double((object_length - reconstruct_box_length) / scanning_step));
inline constexpr std::size_t measure_num = 15;
inline constexpr std::size_t scanning_object_pixel_num = reconstruct_box_pixel_num * (double)(object_length / reconstruct_box_length);

inline constexpr Length distance = 100.0_mm;  // Distance between the object to the detector

inline constexpr Length lambda = 0.5_nm;
inline constexpr WaveNumber k = 2.0 * M_PI / lambda;

inline constexpr PhotonFluxDensity epsilon = 1.0e-20 * dens_unit;
inline constexpr double alpha = 1.0;
inline constexpr double beta = 1.0;

}  // namespace Constant
