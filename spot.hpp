#ifndef SPOT_HPP
#define SPOT_HPP

#include <vector>

class Star;

class Spot {
public:
    Spot(Star* star, double latitude, double longitude, double fillfactor, bool plage, bool mortal);
    std::vector<double> get_ccf(const double time) const;
    double get_flux(const double time) const;
    bool alive(const double time) const;
    bool collides_with(const Spot& other) const;
    double intensity, temperature;

//private:
    Star* star;
    double latitude, longitude, radius, max_radius;
    double time_appear, time_disappear;
    bool plage, mortal;
};

#endif