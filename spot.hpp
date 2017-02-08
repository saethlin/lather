#ifndef SPOT_HPP
#define SPOT_HPP
#include <vector>
#include <math.h>
#include "star.hpp"
#include "profile.hpp"



class Point {
public:
    Point(double x, double y, double z) :
            x(x), y(y), z(z) {}

    Point rotate_x(double angle) {
        double new_y = y*cos(angle) + z*sin(angle);
        double new_z = y*sin(angle) + z*cos(angle);
        return Point(x, new_y, new_z);
    }

    Point rotate_y(double angle) {
        double new_z = z*cos(angle) + x*sin(angle);
        double new_x = z*sin(angle) + x*cos(angle);
        return Point(new_x, y, new_z);
    }

    Point rotate_z(double angle) {
        double new_x = x*cos(angle) + y*sin(angle);
        double new_y = x*sin(angle) + y*cos(angle);
        return Point(new_x, new_y, z);
    }

    double x, y, z;
};


inline double clamp(double lower, double num, double upper) {
    return num <= lower ? lower : num >= upper ? upper : num;
}



class Spot {
public:
    Spot(Star star, double latitude, double longitude, double fillfactor, bool plage, size_t spotResolution);
    bool isVisible(double phase);
    bool isVisible2(double phase);
    void scan(double phase, double& flux, std::vector<double>& profile, double wavelength, bool observeRV);
    double intensity, temperature;

private:
    double latitude, longitude;
    Star star;
    double size;
    bool plage;
    double matrixSpot[3][3];
    std::vector<point> initialCoordinates;
    int iminy, imaxy, iminz, imaxz;
};
#endif