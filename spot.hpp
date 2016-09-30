#ifndef SPOT_HPP
#define SPOT_HPP
#include <vector>
#include <stdbool.h>
#include "star.hpp"


class Spot {
    public:
        Spot(Star star, double latitude, double longitude, double size, bool plage, int spotResolution);

        bool isVisible(double phase);

        void scan(double phase, double& flux, std::vector<double>& profile, bool observeRV);

        Star star;
        double latitude;
        double longitude;
        double size;
        unsigned int spotResolution;
        bool plage;
        double spotTemp, intensity;
        double matrixSpot[3][3];
        std::vector<std::vector<double> > initialCoordinates;
        unsigned int iminy, imaxy, iminz, imaxz;
};
#endif