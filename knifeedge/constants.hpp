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

inline constexpr Length lambda = 0.5_nm;
inline constexpr WaveNumber k = 2.0 * M_PI / lambda;

inline constexpr std::size_t image_pixel_num = 400;
inline constexpr std::size_t detector_pixel_num = 201;

inline constexpr Length detector_pixel_size = 50.0_um;
inline constexpr Length detector_length = detector_pixel_num * detector_pixel_size;

inline constexpr Length knife_step = 50.0_nm;
inline constexpr std::size_t knife_num = 21;
inline constexpr Length probe_length = knife_step * knife_num;
inline constexpr std::size_t probe_num = detector_pixel_num;

inline constexpr double alpha = 0.5;

}  // namespace Constant
