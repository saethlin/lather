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


void normalize(std::vector<double>& vec);


class Star {
public:
    Star() {}
    Star(double radius, double period, double inclination, double temperature, double spotTempDiff,
                    double limbLinear, double limbQuadratic, size_t gridSize);
    double get_limb_integral(const double z_upper, const double z_lower, const double y) const;
    double limb_brightness(const double r_cos) const;
    std::vector<double>& active_profile(const double y);
    std::vector<double>& quiet_profile(const double y);

    friend class Spot;

    double inclination;
    double period;
    double temperature, spotTempDiff;
    double limbLinear, limbQuadratic;
    double intensity;
    double grid_interval;
    double fluxQuiet, zero_rv;
    Profile profileQuiet, profileActive;
    std::vector<double> integrated_ccf;
    std::vector<double> fit_result;

private:
    double equatorial_velocity;

};


#endif