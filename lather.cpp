#include "simulation.hpp"
#include <iostream>


int main() {

    clock_t begin = clock();

    Simulation simulation("/home/ben/lather/config.cfg");

    std::vector<double> time(1000);
    std::iota(time.begin(), time.end(), 0);
    for (auto & val : time) {
        val = val/(double)time.size() * 25.0;
    }

    auto rv = simulation.observe_rv(time, 5000e-10, 5001e-10);


    for (const auto& val : rv) {
        std::cout << val << '\n';
    }

    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //std::cout << "Run took: " << elapsed_secs << std::endl;

    return 0;
}
