#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include "profile.hpp"
#include "fitrv.hpp"
#include "star.hpp"

const double pi = M_PI;
const double solarRadius = 696000.0;

const double c = 299792458;
const double h = 6.62606896e-34;
const double k_b = 1.380e-23;


double planck(double wavelength, double temperature) {
    return 2*h*c*c*1./(wavelength*wavelength*wavelength*wavelength*wavelength*(expm1((h*c)/(wavelength*k_b*temperature))));
}


Star::Star() {}


Star::Star(double radius, double period, double inclination, double temperature, double spotTempDiff,
           double limbLinear, double limbQuadratic, unsigned int gridSize) {
    inclination *= pi/180.0;
    this -> prot = period;
    this -> vrot = (2.0 * pi * solarRadius)/(period * 86400.0);
    this -> inclination = inclination;
    this -> temperature = temperature;
    this -> spotTempDiff = spotTempDiff;
    this -> limbLinear = limbLinear;
    this -> limbQuadratic = limbQuadratic;
    this -> gridSize = gridSize;
    this -> intensity = planck(5293.4115e-10, temperature);

    // Setup for profiles
    std::vector<double> rv;
    ccfQuiet = std::vector<double>(0);
    std::vector<double> ccfActive;

    std::string filename;
    if (vrot < 10) filename = "resources/solarccfhires.txt";
    else filename = "resources/solarccfhires.txt";

    std::ifstream ifs(filename);
    std::string line;
    double num;

    // Skip first two header lines
    getline(ifs, line);
    getline(ifs, line);

    while (getline(ifs, line)) {
        std::istringstream stream(line);
        stream >> num;
        rv.push_back(num);

        stream >> num;
        ccfQuiet.push_back(num);

        stream >> num;
        ccfActive.push_back(num);
    }

    profileQuiet = Profile(rv, ccfQuiet);
    profileActive = Profile(rv, ccfActive);

    unsigned int iy, iz, k;
    double y, z;
    double v_shift, r_cos, limbSum;
    std::vector<double> ccfShifted;
    fit_result = std::vector<double>(4);
    ccfQuiet = std::vector<double>(profileQuiet.size);
    limb = std::vector<double>(gridSize*gridSize);
    double rSquared;

    fluxQuiet = 0;
    for (iy = 0; iy <= gridSize; iy++) {
        y = -1.0+iy*2.0/gridSize;

        v_shift = y * vrot * sin(inclination);
        ccfShifted = profileQuiet.shift(v_shift);

        limbSum = 0.0;

        for (iz = 0; iz <= gridSize; iz++) {
            z = -1.0+iz*2.0/gridSize;

            rSquared = y*y + z*z;

            if (rSquared <= 1) {
                r_cos = sqrt(1 - rSquared);
                limbSum += 1 - limbLinear*(1.-r_cos) - limbQuadratic*(1-r_cos)*(1-r_cos);
                //cachePos = rSquared*gridSize*gridSize;
                //limb[cachePos] = 1 - limbLinear*(1.-r_cos) - limbQuadratic*(1-r_cos)*(1-r_cos);
                //limbSum += limb[cachePos];
            }
        }

        for (k = 0; k < profileQuiet.size; k++) {
            ccfQuiet[k] += ccfShifted[k] * limbSum;
        }

        fluxQuiet += limbSum;
    }

    // Compute the rv that will be fitted with no spots visible.
    std::vector<double> normProfile(ccfQuiet);
    normalize(normProfile);

    fit_result[0] = normProfile[normProfile.size()/2] - normProfile[0];
    fit_result[1] = profileQuiet.rv[normProfile.size()/2];
    fit_result[2] = 2.71; // TODO Remove this magic number
    fit_result[3] = normProfile[0];

    fit_rv(profileQuiet.rv, normProfile, fit_result);
    zero_rv = fit_result[1];
}
