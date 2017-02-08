#include "simulation.hpp"
#include <iostream>


int main() {
    clock_t begin = clock();
    Simulation simulation("/home/ben/lather/config.cfg");
    size_t npoints = 25;
    std::vector<double> time(npoints);
    for (auto i = 0; i < time.size(); i++) {
        time[i] = (double)i/time.size()*25.05;
    }
    std::vector<double> flux(npoints);
    std::vector<double> rv(npoints);
    simulation.observe(time, flux, rv, 5000e-10, true);
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << elapsed_secs << std::endl;

    for (const auto& val : rv) {
        std::cout << val << std::endl;
    }

    return 0;
}
