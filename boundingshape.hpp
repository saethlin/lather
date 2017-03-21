#ifndef LATHER_BOUNDINGSHAPE_H
#define LATHER_BOUNDINGSHAPE_H


#include "point.hpp"
#include "spot.hpp"
#include <cmath>
#include <iostream>


struct bounds {
    bounds(double val1, double val2) {
        lower = std::min(val1, val2);
        upper = std::max(val1, val2);
    }
    double lower, upper;
};


class Spot;


class BoundingShape {
public:
    BoundingShape(const Spot& spot, const double phase);
    bounds y_bounds() const;
    bounds z_bounds(const double y) const;
    bounds z_bounds_edge(const double y) const;
    bounds z_bounds_edge2(const double y) const;

private:
    bool on_spot(const double y, const double z) const;
    Point center, circle_center, a, b;
    double radius, grid_interval;
    bool visible, is_on_edge;
};


#endif
