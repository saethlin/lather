#include "simulation.hpp"
#include <iostream>


int main() {

    clock_t begin = clock();

    Simulation simulation("/home/ben/lather/config.cfg");

    std::vector<double> time(25);
    std::iota(time.begin(), time.end(), 0);

    auto rv = simulation.observe_flux(time, 4000e-10, 5000e-10);

    for (const auto& val : rv) {
        std::cout << val << '\n';
    }

    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "Run took: " << elapsed_secs << std::endl;

    return 0;
}
