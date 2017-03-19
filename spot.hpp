#ifndef SPOT_HPP
#define SPOT_HPP


#include "star.hpp"
#include "profile.hpp"
#include "point.hpp"
#include "boundingshape.hpp"
#include <vector>
#include <math.h>
#include <numeric>


class Spot {
public:
    Spot(Star* star, double latitude, double longitude, double fillfactor, bool plage);
    std::vector<double> get_ccf(const double time) const;
    double get_flux(const double time) const;
    bool alive(const double time) const;
    double intensity, temperature;

//private:
    Star* star;
    double latitude, longitude, radius;
    double time_appear, time_disappear;
    bool plage;
};

#endif