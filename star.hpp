#ifndef STAR_HPP
#define STAR_HPP
#include <vector>
#include "profile.hpp"

void normalize(std::vector<double>& vec);

class Star {
    public:
        Star();
        Star(double radius, double period, double inclination, double temperature, double spotTempDiff,
                        double limbLinear, double limbQuadratic, unsigned int gridSize);

        double vrot;
        double prot;
        double inclination;
        double temperature, spotTempDiff;
        double limbLinear, limbQuadratic;
        double intensity;
        unsigned int gridSize;
        int spotResolution;
        Profile profileQuiet;
        Profile profileActive;
        double fluxQuiet;
        std::vector<double> ccfQuiet;
        std::vector<double> limb;
        double zero_rv;
        std::vector<double> fit_result;
};

double planck(double wavelength, double temperature);

#endif