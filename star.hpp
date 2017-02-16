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
    Star();
    Star(double radius, double period, double inclination, double temperature, double spotTempDiff,
                    double limbLinear, double limbQuadratic, size_t gridSize);
    double limb_integral(double, double, const double) const;
    double limb_brightness(const double r_cos) const;
    std::vector<double>& active_profile(const double y);
    std::vector<double>& quiet_profile(const double y);
    bool load_cache();
    void save_cache();

    double inclination;
    double period;
    double temperature, spotTempDiff;
    double limbLinear, limbQuadratic;
    double intensity;
    double grid_interval;
    Profile profileQuiet;
    Profile profileActive;
    double fluxQuiet;
    double zero_rv;
    std::vector<double> integrated_ccf;
    std::vector<double> fit_result;

private:
    double equatorial_velocity;

};

double planck(double wavelength, double temperature);

#endif