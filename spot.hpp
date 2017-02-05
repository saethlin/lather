#ifndef SPOT_HPP
#define SPOT_HPP
#include <vector>
#include <math.h>
#include "star.hpp"
#include "profile.hpp"



class Spot {
    public:
        Spot(Star star, double latitude, double longitude, double size, bool plage, int spotResolution);

        bool isVisible(double phase);

        void scan(double phase, double& flux, std::vector<double>& profile, double wavelength, bool observeRV);

        Star star;
        double latitude;
        double longitude;
        double size;
        int spotResolution;
        bool plage;
        double spotTemp, intensity;
        double matrixSpot[3][3];
        std::vector<std::vector<double> > initialCoordinates;
        int iminy, imaxy, iminz, imaxz;
};
#endif