#pragma once

#include <unit/double.hpp>

inline namespace Constant
{


using namespace Unit;

using WaveNumber = decltype(1.0 / Length{});
using PhotonFluxDensity = DimensionLessType;
using WaveAmplitude = decltype(PhotonFluxDensity{}.sqrt());

inline const WaveAmplitude amp_unit{1.0};
inline const PhotonFluxDensity dens_unit{1.0};

inline constexpr std::size_t image_pixel_num = 847;
inline constexpr std::size_t detector_pixel_num = 201;
inline constexpr std::size_t reconstruct_box_pixel_num = detector_pixel_num;

inline constexpr Length detector_pixel_size = 50.0_um;
inline constexpr Length detector_length = double(detector_pixel_num) * detector_pixel_size;
inline constexpr Length object_length = 5.5_um;

inline constexpr Length light_stddev = 200.0_nm;
inline constexpr Length reconstruct_box_length = 500.0_nm;
inline constexpr Length reconstruct_box_pixel_size = reconstruct_box_length / (double)reconstruct_box_pixel_num;

inline constexpr Length scanning_step = 100.0_nm;
inline constexpr Length scanning_half_length = 0.5 * (object_length - reconstruct_box_length);
inline constexpr double overlap_ratio = 1.0 - (double)(scanning_step / reconstruct_box_length);

inline constexpr std::size_t measure_num = double((object_length - reconstruct_box_length) / scanning_step) + 1;
inline constexpr std::size_t scanning_object_pixel_num = reconstruct_box_pixel_num * double(object_length / reconstruct_box_length);

inline constexpr Length distance = 100.0_mm;  // Distance between the object to the detector

inline constexpr Length lambda = 0.5_nm;
inline constexpr WaveNumber k = 2.0 * M_PI / lambda;

inline constexpr double alpha = 0.5;
inline constexpr double beta = 0.5;

}  // namespace Constant
