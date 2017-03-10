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
    Simulation() {}
    Simulation(const char* filename);
    Simulation(size_t gridSize);
    void set_star(double radius, double period, double inclination, double temperature, double spot_temp_diff,
                  double linear_limb, double quadratic_limb);
    void add_spot(double latitude, double longitude, double size, bool plage);
    void clear_spots();
    std::vector<double> observe_rv(const std::vector<double>& time, const double wavelength_min, const double wavelength_max);
    std::vector<double> observe_flux(const std::vector<double>& time, const double wavelength_min, const double wavelength_max);

private:
    size_t grid_size;
    Star star;
    std::vector<Spot> spots;
};
#endif
