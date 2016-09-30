#include <iostream>
#include <string>
#include <vector>
#include <stdbool.h>
#include "inih/cpp/INIReader.h"
#include "star.hpp"
#include "spot.hpp"
#include "fitrv.hpp"
#include "fitsim.hpp"
#include "simulation.hpp"

const double pi = M_PI;


void normalize(std::vector<double>& vec) {
    unsigned int i;
    double max = vec[0];
    for (i = 0; i < vec.size(); i++) {
        if (vec[i] > max) max = vec[i];
    }
    for (i = 0; i < vec.size(); i++) {
        vec[i] = vec[i]/max;
    }
}


Simulation::Simulation(unsigned int gridSize, unsigned int spotResolution) {
    this -> gridSize = gridSize;
    this -> spotResolution = spotResolution;
}


Simulation::Simulation(const char* filename) {

    double radius, prot, inclination, temperature, spot_temp_diff, limbLinear, limbQuadratic;
    double latitude, longitude, size;
    bool plage;

    INIReader reader(filename);
    if (reader.ParseError() < 0) {
        throw std::exception();
    }

    this -> gridSize = reader.GetInteger("simulation", "grid_resolution", 100);
    this -> spotResolution = reader.GetInteger("simulation", "spot_resolution", 100);

    radius = reader.GetReal("star", "radius", 1.0);
    prot = reader.GetReal("star", "prot", 25.05);
    inclination = reader.GetReal("star", "inclination", 90.0);
    temperature = reader.GetReal("star", "Tstar", 5778.0);
    spot_temp_diff = reader.GetReal("star", "Tdiff_spot", 663.0);
    limbLinear = reader.GetReal("star", "limb1", 0.29);
    limbQuadratic = reader.GetReal("star", "limb2", 0.34);

    setStar(radius, prot, inclination, temperature, spot_temp_diff, limbLinear, limbQuadratic);

    std::set<std::string> sections = reader.GetSections();
    std::set<std::string>::iterator sectionName = sections.begin();
    for (sectionName = sections.begin(); sectionName != sections.end(); ++sectionName) {

        std::string name = *sectionName;
        if ((*sectionName).substr(0, 4) == "spot") {

            latitude = reader.GetReal(name, "latitude", 0.0);
            longitude = reader.GetReal(name, "longitude", 180.0);
            size = reader.GetReal(name, "size", 0.1);
            plage = reader.GetBoolean(name, "plage", false);

            addSpot(latitude, longitude, size, plage);
        }
    }
}


void Simulation::setStar(double radius, double period, double inclination, double temperature, double spot_temp_diff, double linear_limb, double quadratic_limb) {
    star = Star(radius, period, inclination, temperature, spot_temp_diff, linear_limb, quadratic_limb, gridSize);
}


void Simulation::addSpot(double latitude, double longitude, double size, bool plage) {
    Spot new_spot(star, latitude, longitude, size, plage, spotResolution);
    spots.push_back(new_spot);
}


void Simulation::observe(std::vector<double>& time, std::vector<double>& flux, std::vector<double>& rv, bool observeRV) {
    unsigned int s, t, i;
    double phase;
    double spotFlux = 0;
    bool anyVisible = false;

    std::vector<double> spotProfile(star.profileQuiet.ccf.size());
    std::vector<double> fit_result = star.fit_result;

    for (t = 0; t < time.size(); t++) {
        phase = fmod(time[t], star.prot)/star.prot * 2*pi;

        for (s = 0; s < spots.size(); s++) {
            if (spots[s].isVisible(phase)) {
                anyVisible = true;
                spots[s].scan(phase, spotFlux, spotProfile, observeRV);
            }
        }

        if (anyVisible) {
            anyVisible = false;

            flux[t] = (star.fluxQuiet - spotFlux) / star.fluxQuiet;
            spotFlux = 0.0;

            if (observeRV) {
                // Compute the observed profile and fit the rv: the star's quiet profile minus the spot flux
                for (i = 0; i < spotProfile.size(); i++) {
                    spotProfile[i] = star.ccfQuiet[i] - spotProfile[i];
                }
                normalize(spotProfile);
                fit_rv(star.profileQuiet.rv, spotProfile, fit_result);
                rv[t] = fit_result[1] - star.zero_rv;

                for (i = 0; i < spotProfile.size(); i++) {
                    spotProfile[i] = 0.0;
                }
            }
        }
        else {
            rv[t] = 0.0;
            flux[t] = 1.0;
        }
    }
    normalize(flux);
}


void Simulation::fit(std::vector<double>& time, std::vector<double>& flux) {
    fit_sim(this, time, flux);
}
