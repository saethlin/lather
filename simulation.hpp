#ifndef SIMULATION_HPP
#define SIMULATION_HPP


#include "spot.hpp"
#include "star.hpp"
#include "fitrv.hpp"
#include "inih/cpp/INIReader.h"
#include <vector>


class Simulation {
    public:
        Simulation() {}

        Simulation(const char* filename);

        Simulation(unsigned int gridSize, unsigned int spotResolution);

        void setStar(double radius, double period, double inclination, double temperature, double spot_temp_diff, double linear_limb, double quadratic_limb);

        void addSpot(double latitude, double longitude, double size, bool plage);

        void observe(std::vector<double>& time, std::vector<double>& flux, std::vector<double>& rv, double wavelength, bool observeRV);

        //void fit(std::vector<double>& time, std::vector<double>& flux);

        unsigned int gridSize;
        unsigned int spotResolution;
        Star star;
        std::vector<Spot> spots;
};
#endif
