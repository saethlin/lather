#ifndef SIMULATION_HPP
#define SIMULATION_HPP


#include "spot.hpp"
#include "star.hpp"
#include "fitrv.hpp"
#include "planck.hpp"
#include "inih/cpp/INIReader.h"
#include <vector>
#include <algorithm>


class Simulation {
public:
    //Simulation() {}
    Simulation(const char* filename);
    Simulation(size_t gridSize);
    void setStar(double radius, double period, double inclination, double temperature, double spot_temp_diff, double linear_limb, double quadratic_limb);
    void addSpot(double latitude, double longitude, double size, bool plage);
    void clear_spots();
    std::vector<double> observe_rv(std::vector<double>& time, double wavelength_min, double wavelength_max);
    std::vector<double> observe_flux(std::vector<double>& time, double wavelength_min, double wavelength_max);
    //void fit(std::vector<double>& time, std::vector<double>& flux);

private:
    size_t gridSize;
    Star star;
    std::vector<Spot> spots;
};
#endif
