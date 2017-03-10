#include "point.hpp"
#include <cmath>


void Point::rotate_x(double angle) {
    double new_y = y*cos(angle) - z*sin(angle);
    double new_z = y*sin(angle) + z*cos(angle);
    y = new_y;
    z = new_z;
}


void Point::rotate_y(double angle) {
    double new_z = z*cos(angle) - x*sin(angle);
    double new_x = z*sin(angle) + x*cos(angle);
    z = new_z;
    x = new_x;
}


void Point::rotate_z(double angle) {
    double new_x = x*cos(angle) - y*sin(angle);
    double new_y = x*sin(angle) + y*cos(angle);
    x = new_x;
    y = new_y;
}
