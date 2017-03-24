#ifndef SIMULATION_HPP
#define SIMULATION_HPP


#include "spot.hpp"
#include "star.hpp"
#include "fitrv.hpp"
#include "planck.hpp"
#include "inih/cpp/INIReader.h"
#include <vector>
#include <algorithm>
#include <Magick++.h>
#include <sstream>
#include <iomanip>


class Simulation {
public:
    Simulation() {}
    Simulation(const char* filename);
    void set_star(const double radius, const double period, const double inclination, const double temperature,
                  const double spot_temp_diff, const double linear_limb, const double quadratic_limb, const size_t grid_size);
    void add_spot(const double latitude, const double longitude, const double size, const bool plage, const bool mortal);
    void clear_spots();
    std::vector<double> observe_rv(const std::vector<double>& time, const double wavelength_min, const double wavelength_max);
    std::vector<double> observe_flux(const std::vector<double>& time, const double wavelength_min, const double wavelength_max);
    void draw(const double time, const int i) const;
    std::vector<uint8_t> draw_rgba(const double time);

    Simulation(Simulation const&) = delete;
    Simulation& operator=(Simulation const&) = delete;


private:
    void check_fill_factor(double time);
    Star star;
    std::vector<Spot> spots;
    double dynamic_fill_factor;
};
#endif
