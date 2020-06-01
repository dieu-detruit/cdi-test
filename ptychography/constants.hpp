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
inline constexpr std::size_t object_pixel_num = image_pixel_num;
inline constexpr std::size_t detector_pixel_num = 201;
inline constexpr std::size_t light_pixel_num = detector_pixel_num;

inline constexpr Length detector_pixel_size = 50.0_um;
inline constexpr Length detector_length = double(detector_pixel_num) * detector_pixel_size;
inline constexpr Length object_length = 6.0_um;

inline constexpr Length object_pixel_size = object_length / double(object_pixel_num);
inline constexpr Length light_stddev = 150.0_nm;
inline constexpr Length light_length = 2.0 * light_stddev;

inline constexpr Length scanning_step = 100.0_nm;
inline constexpr std::size_t image_num = std::size_t(double((object_length - light_length) / scanning_step));
inline constexpr Length scanning_half_length = scanning_step * double((image_num - 1) / 2);
inline constexpr double overlap_ratio = 1.0 - (double)(scanning_step / light_length);

inline constexpr Length distance = 100.0_mm;  // Distance between the object to the detector

inline constexpr Length lambda = 0.5_nm;
inline constexpr WaveNumber k = 2.0 * M_PI / lambda;

inline constexpr double alpha = 0.85;
inline constexpr double beta = 0.85;

}  // namespace Constant
