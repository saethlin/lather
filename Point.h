#ifndef LATHER_POINT_H
#define LATHER_POINT_H


// y and z max (the spot's radius) are sqrt(2*x-x*x)

struct point {
    double x, y, z;
};


class PointVector {
public:
    void push_back(point);
    void rotate();
    PointVector get_rotated(angle);
};


void PointVector::rotate_x(rotmatrix m) {

}


#endif