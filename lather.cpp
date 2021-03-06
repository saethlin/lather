#include "simulation.hpp"
#include <iostream>


int main() {
    clock_t begin = clock();

    Simulation simulation("/home/ben/research/lather_experiments/sun.cfg");

    std::vector<double> time(100);
    std::iota(time.begin(), time.end(), 0);

    for (auto& val : time) {
        val = val/(double)time.size() * 25.0;
    }

    auto rv = simulation.observe_rv(time, 4000e-10, 5000e-10);

    for (const auto& val : rv) {
        //std::cout << val.rv << '\n';
        //std::cout << val.bisector.size() << '\n';
    }

    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << "Run took: " << elapsed_secs << std::endl;

    return 0;
}
