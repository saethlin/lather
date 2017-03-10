#ifndef LATHER_POINT_H
#define LATHER_POINT_H

class Point {
public:
    Point() {};
    Point(double x, double y, double z) :
            x(x), y(y), z(z) {}
    void rotate_x(double angle);
    void rotate_y(double angle);
    void rotate_z(double angle);

    double x, y, z;
};


#endif