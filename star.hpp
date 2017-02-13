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
        double limb_integral(double, double, double);

        double vrot;
        double period;
        double inclination;
        double temperature, spotTempDiff;
        double limbLinear, limbQuadratic;
        double intensity;
        double grid_interval;
        Profile profileQuiet;
        Profile profileActive;
        double fluxQuiet = 0;
        std::vector<double> integrated_ccf;
        std::vector<double> limb;
        double zero_rv;
        std::vector<double> fit_result;
        bool analytic;
};

double planck(double wavelength, double temperature);

#endif