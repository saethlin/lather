#ifndef STAR_HPP
#define STAR_HPP

#include <vector>
#include <random>
#include "profile.hpp"

void normalize(std::vector<double>& vec);


class Star {
public:
    Star() {}
    Star(const double radius, const double period, const double inclination, const double temperature,
         const double spot_temp_diff, const double limb_linear, const double limb_quadratic, const size_t gridSize);
    double get_limb_integral(const double z_upper, const double z_lower, const double y) const;
    double limb_brightness(const double r_cos) const;
    std::vector<double>& active_profile(const double y) const;
    std::vector<double>& quiet_profile(const double y) const;

    double period, inclination, temperature, spot_temp_diff, limb_linear, limb_quadratic, grid_interval;
    double intensity, flux_quiet, zero_rv;
    Profile profile_quiet, profile_active;
    std::vector<double> integrated_ccf;
    std::vector<double> fit_result;
    std::vector<float> image;

    double equatorial_velocity;
    std::mt19937 generator;

private:

};


#endif