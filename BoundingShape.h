#ifndef LATHER_BOUNDINGSHAPE_H
#define LATHER_BOUNDINGSHAPE_H

#include "spot.hpp"


class Spot;


struct bounds {
    double lower, upper;
};


class BoundingShape {
public:
    BoundingShape(const Spot& spot, const double phase);
    bool is_visible() const;
    bounds get_y_bounds() const;
    bounds get_z_bounds(const double y) const;

private:
    Point center, circle_center, a, b;
    double size, radius, grid_interval;
};


#endif
