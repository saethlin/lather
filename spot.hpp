#ifndef SPOT_HPP
#define SPOT_HPP
#include <vector>
#include <math.h>
#include "star.hpp"
#include "profile.hpp"


class Point {
public:
    Point() {};
    Point(double x, double y, double z) :
            x(x), y(y), z(z) {}

    Point rotate_x(double angle) {
        double new_y = y*cos(angle) - z*sin(angle);
        double new_z = y*sin(angle) + z*cos(angle);
        return Point(x, new_y, new_z);
    }

    Point rotate_y(double angle) {
        double new_z = z*cos(angle) - x*sin(angle);
        double new_x = z*sin(angle) + x*cos(angle);
        return Point(new_x, y, new_z);
    }

    Point rotate_z(double angle) {
        double new_x = x*cos(angle) - y*sin(angle);
        double new_y = x*sin(angle) + y*cos(angle);
        return Point(new_x, new_y, z);
    }

    Point rotate_axis(double angle, Point axis) {
        double u = axis.x;
        double v = axis.y;
        double w = axis.z;

        double new_x = u*(u*x + v*y + w*z)*(1-cos(angle)) + x*cos(angle) + (-w*y + v*z)*sin(angle);
        double new_y = v*(u*x + v*y + w*z)*(1-cos(angle)) + y*cos(angle) + (w*x - u*z)*sin(angle);
        double new_z = w*(u*x + v*y + w*z)*(1-cos(angle)) + z*cos(angle) + (-v*x + u*y)*sin(angle);
        return Point(new_x, new_y, new_z);
    }

    double x=0, y=0, z=0;
};


inline double clamp(double lower, double num, double upper) {
    return num <= lower ? lower : num >= upper ? upper : num;
}


struct spot_bounds {
    bool visible;
    double miny, maxy, minz, maxz;
};


class Spot {
public:
    Spot(Star star, double latitude, double longitude, double fillfactor, bool plage, size_t spotResolution);
    spot_bounds get_bounds_at(double phase);
    bool isVisible2(double phase);
    bool is_on_spot(double x, double y, double z);
    std::vector<double> get_ccf(double phase, double wavelength);
    double get_flux(double phase, double wavelength);
    double intensity, temperature;

private:
    Star star;
    double opening_angle, radius;
    double radius_squared;
    double latitude, longitude;
    double size;
    bool plage;
    double matrixSpot[3][3];
    std::vector<Point> initialCoordinates;
    int iminy, imaxy, iminz, imaxz;
    double miny, maxy, minz, maxz;
    Point center;
};
#endif