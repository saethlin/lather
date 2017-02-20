#ifndef SPOT_HPP
#define SPOT_HPP
#include <vector>
#include <math.h>
#include "star.hpp"
#include "profile.hpp"
#include "Point.h"
#include "BoundingShape.h"


struct spot_bounds {
    bool visible;
    double miny, maxy, minz, maxz;
};


class Spot {
public:
    Spot(Star* star, double latitude, double longitude, double fillfactor, bool plage);
    spot_bounds get_bounds_at(double phase);
    bool is_on_spot(double x, double y, double z);
    std::vector<double> get_ccf(double phase);
    std::vector<double> get_ccf_dev(double phase);
    double get_flux(double phase);
    double intensity, temperature;

//private:
    Star* star;
    double latitude, longitude, size, radius_squared;
    bool plage;
    Point center;
};

#endif