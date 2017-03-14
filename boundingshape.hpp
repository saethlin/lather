#ifndef LATHER_BOUNDINGSHAPE_H
#define LATHER_BOUNDINGSHAPE_H


#include "point.hpp"
#include "spot.hpp"


struct bounds {
    double lower, upper;
};


class Spot;


class BoundingShape {
public:
    BoundingShape(const Spot& spot, const double phase);
    bool is_visible() const;
    bounds y_bounds() const;
    bounds z_bounds(const double y) const;

private:
    bool on_spot(const double y, const double z) const;
    Point center, circle_center, a, b;
    double size, radius, grid_interval, time;
};


#endif
