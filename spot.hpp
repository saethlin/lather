#ifndef SPOT_HPP
#define SPOT_HPP


#include "star.hpp"
#include "profile.hpp"
#include "Point.h"
#include "BoundingShape.h"
#include <vector>
#include <math.h>


class Spot {
public:
    Spot(Star* star, double latitude, double longitude, double fillfactor, bool plage);
    std::vector<double> get_ccf(double phase);
    double get_flux(double phase);
    double intensity, temperature;

//private:
    Star* star;
    double latitude, longitude, size;
    bool plage;
};

#endif