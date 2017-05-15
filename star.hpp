#ifndef STAR_HPP
#define STAR_HPP


#include "profile.hpp"
#include "fitrv.hpp"
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include <functional>
#include <gsl/gsl_integration.h>
#include <Magick++.h>


void normalize(std::vector<double>& vec);


class Star {
public:
    Star() {}
    Star(const double radius, const double period, const double inclination, const double temperature,
         const double spot_temp_diff, const double limb_linear, const double limb_quadratic, const size_t gridSize);
    double get_limb_integral(const double z_upper, const double z_lower, const double y) const;
    double limb_brightness(const double r_cos) const;
    double get_velocity(const double y, const double z) const;
    double Star::limb_path_func(const double v, void* args) const;

    double period, inclination, temperature, spot_temp_diff, limb_linear, limb_quadratic, grid_interval, velocity_interval;
    double intensity, flux_quiet, zero_rv;
    Profile profile_quiet, profile_active;
    std::vector<double> integrated_ccf;
    std::vector<double> fit_result;
    std::vector<float> image;

    double equatorial_velocity;
    double diff_a, diff_b, diff_c;
};


#endif