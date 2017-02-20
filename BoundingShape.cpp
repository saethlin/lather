#include "BoundingShape.h"


BoundingShape::BoundingShape(const Spot* spot, const double phase) {
    this->spot = spot;

    double theta = phase + spot->longitude;
    double phi = M_PI_2 - spot->latitude;
    center = Point(sin(phi)*cos(theta),
                   sin(phi)*sin(theta),
                   cos(phi));
    center.rotate_y(spot->star->inclination - M_PI_2);

    double depth = sqrt(1-spot->size*spot->size);
    circle_center = {center.x*depth, center.y*depth, center.z*depth};

    auto h = 1-depth;
    radius = sqrt(2*h - h*h);

    double a_x = 0;
    double a_y = -circle_center.z/sqrt(circle_center.y*circle_center.y + circle_center.z*circle_center.z);
    double a_z = sqrt(1-a_y*a_y);
    a = {a_x, a_y, a_z};

    double b_x = circle_center.y*a_z - circle_center.z*a_y;
    double b_y = circle_center.z*a_x - circle_center.x*a_z;
    double b_z = circle_center.x*a_y - circle_center.y*a_x;

    b = {b_x, b_y, b_z};
}


bool BoundingShape::is_visible() const {
    return center.x > -sqrt(2*spot->size);
}


bounds BoundingShape::get_y_bounds() const {
    double theta_y_max = -2 * atan(a.y/b.y - (sqrt(a.y*a.y + b.y*b.y)/b.y));
    double theta_y_min = -2 * atan(a.y/b.y + (sqrt(a.y*a.y + b.y*b.y)/b.y));

    double y_max = circle_center.y + radius*(cos(theta_y_max)*a.y + sin(theta_y_max)*b.y);
    double y_min = circle_center.y + radius*(cos(theta_y_min)*a.y + sin(theta_y_min)*b.y);

    y_min = floor(y_min/spot->star->grid_interval)*spot->star->grid_interval;
    y_max = ceil(y_max/spot->star->grid_interval)*spot->star->grid_interval;

    return {y_min, y_max};
}


bounds BoundingShape::get_z_bounds(const double y) const {

    double interior = a.y*a.y*radius*radius + b.y*b.y*radius*radius - circle_center.y*circle_center.y + 2.*circle_center.y*y - y*y;
    interior = std::max(interior, 0.0);

    double theta_z_max = 2. * atan((-b.y*radius + sqrt(interior))/
                                  (-a.y*radius + center.y - y));
    double theta_z_min = -2. * atan((b.y*radius + sqrt(interior))/
                                  (-a.y*radius + center.y - y));

    double z_max = circle_center.z + radius*(cos(theta_z_max)*a.z + sin(theta_z_max)*b.z);
    double z_min = circle_center.z + radius*(cos(theta_z_min)*a.z + sin(theta_z_min)*b.z);

    z_min = floor(z_min/spot->star->grid_interval)*spot->star->grid_interval;
    z_max = ceil(z_max/spot->star->grid_interval)*spot->star->grid_interval;

    return {z_min, z_max};
}
