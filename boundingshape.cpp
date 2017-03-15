#include "boundingshape.hpp"


BoundingShape::BoundingShape(const Spot& spot, const double time) {
    grid_interval = spot.star->grid_interval;
    size = spot.size;
    this-> time = time;

    const double phase = fmod(time, spot.star->period) / spot.star->period * 2 * M_PI;
    const double theta = phase + spot.longitude;
    const double phi = M_PI_2 - spot.latitude;
    center = Point(sin(phi)*cos(theta),
                   sin(phi)*sin(theta),
                   cos(phi));
    center.rotate_y(spot.star->inclination - M_PI_2);

    const double depth = sqrt(1-size*size);
    circle_center = {center.x*depth, center.y*depth, center.z*depth};

    const double h = 1-depth;
    radius = sqrt(2*h - h*h);

    const double a_x = 0;
    const double a_y = -circle_center.z/sqrt(circle_center.y*circle_center.y + circle_center.z*circle_center.z);
    const double a_z = sqrt(1-a_y*a_y);
    a = {a_x, a_y, a_z};

    const double b_x = circle_center.y*a_z - circle_center.z*a_y;
    const double b_y = circle_center.z*a_x - circle_center.x*a_z;
    const double b_z = circle_center.x*a_y - circle_center.y*a_x;

    b = {b_x, b_y, b_z};

    visible = center.x > -radius;

    // Look for the smallest x value that the boundary touches to determine if the spot is entirely on the visible hemisphere
    double theta_x_max = -2 * atan((a.x - sqrt(a.x * a.x + b.x * b.x)) / b.x);
    double theta_x_min = -2 * atan((a.x + sqrt(a.x * a.x + b.x * b.x)) / b.x);

    double x_1 = circle_center.x + radius*(cos(theta_x_max)*a.x + sin(theta_x_max)*b.x);
    double x_2 = circle_center.x + radius*(cos(theta_x_min)*a.x + sin(theta_x_min)*b.x);
    
    double x_min = std::min(x_1, x_2);
    double x_max = std::max(x_1, x_2);

    is_on_edge = x_min < 0 || x_max < 0; // Check if signs are different
}


bool BoundingShape::is_visible() const {
    return visible;
}


bounds BoundingShape::y_bounds() const {

    if (not visible) {
        return {0.0, 0.0};
    }

    double theta_y_max = -2 * atan((a.y - sqrt(a.y * a.y + b.y * b.y)) / b.y);
    double theta_y_min = -2 * atan((a.y + sqrt(a.y * a.y + b.y * b.y)) / b.y);
    if (b.y == 0) {
        theta_y_max = M_PI;
        theta_y_min = 0;
    }

    double y_max = circle_center.y + radius*(cos(theta_y_max)*a.y + sin(theta_y_max)*b.y);
    double y_min = circle_center.y + radius*(cos(theta_y_min)*a.y + sin(theta_y_min)*b.y);

    y_max = ceil(y_max/grid_interval)*grid_interval;
    y_min = floor(y_min/grid_interval)*grid_interval;

    // Check if the spot is past the visible edge of the star
    const double x_max = circle_center.x + radius*(cos(theta_y_max)*a.y + sin(theta_y_max)*b.y);
    const double x_min = circle_center.x + radius*(cos(theta_y_min)*a.y + sin(theta_y_min)*b.y);

    if (x_min < 0 && x_max < 0) {
        return {0.0, 0.0};
    }

    // Spot wraps around the right edge of the star
    if (x_max < 0.0) {
        y_max = 1.0;
    }

    // Spot wraps around the left edge of the star
    if (x_min < 0.0) {
        y_min = -1.0;
    }

    return {y_min, y_max};
}


bounds BoundingShape::z_bounds(const double y) const {

    if (is_on_edge) {
        return z_bounds_edge(y);
    }

    double interior = a.y*a.y*radius*radius + b.y*b.y*radius*radius - circle_center.y*circle_center.y + 2.*circle_center.y*y - y*y;
    interior = std::max(interior, 0.0);

    const double theta_z_max = 2. * atan((-b.y*radius + sqrt(interior))/
                                  (-a.y*radius + center.y - y));
    const double theta_z_min = -2. * atan((b.y*radius + sqrt(interior))/
                                  (-a.y*radius + center.y - y));

    double z_max = circle_center.z + radius*(cos(theta_z_max)*a.z + sin(theta_z_max)*b.z);
    double z_min = circle_center.z + radius*(cos(theta_z_min)*a.z + sin(theta_z_min)*b.z);

    return {z_min, z_max};
}


bool on_star(const double y, const double z) {
    return (y*y + z*z) <= 1;
}


bool BoundingShape::on_spot(const double y, const double z) const {
    if (!on_star(y, z)) {
        return false;
    }
    double x = sqrt(1-(y*y + z*z));
    double distance_squared = (y-center.y)*(y-center.y) + (z-center.z)*(z-center.z) + (x-center.x)*(x-center.x);
    return (distance_squared <= (radius*radius));
}


bounds BoundingShape::z_bounds_edge(const double y) const {

    // Set both bounds to invalid values so we can detect if a spot was found
    double z_max = 2.0;
    double z_min = 2.0;

    for (double z = center.z+radius; z > center.z-radius; z-=grid_interval) {
        if (on_spot(y, z)) {
            z_max = z;
            break;
        }
    }

    for (double z = center.z-radius; z < center.z+radius; z+=grid_interval) {
        if (on_spot(y, z)) {
            z_min = z;
            break;
        }
    }

    if (z_min > 1.0 || z_max > 1.0) {
        z_min = 0.0;
        z_max = 0.0;
    }

    return {z_min, z_max};
}

/*0.470501 0
0.470127 0.470125
0.469751 0.469749
*/