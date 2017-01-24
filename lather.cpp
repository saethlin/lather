#include <string>
#include <sstream>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include "simulation.hpp"


int main() {

    Simulation simulation("/home/ben/lather/config.cfg");
    std::vector<double> time(1000);
    for (unsigned int i = 0; i < time.size(); i++) {
        time[i] = (float)i/time.size()*25.05;
    }
    std::vector<double> flux(1000);
    std::vector<double> rv(1000);
    std::cout << "data made " << std::endl;
    simulation.observe(time, flux, rv, true);

    return 0;
}
