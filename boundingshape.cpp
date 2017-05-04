#include "boundingshape.hpp"


BoundingShape::BoundingShape(const Spot& spot, const double time) {
    grid_interval = spot.star->grid_interval;
    radius = spot.radius;
    max_radius = radius;

    if (spot.mortal) {
        const double lifetime = spot.time_disappear - spot.time_appear;
        const double growth_time = 0.1 * lifetime;

        if (fabs(time - spot.time_appear) < growth_time) {
            radius *= fabs(time - spot.time_appear) / growth_time;
        } else if (fabs(time - spot.time_disappear) < growth_time) {
            radius *= fabs(time - spot.time_disappear) / growth_time;
        }
    }

    const double phase = fmod(time, spot.star->period) / spot.star->period * 2 * M_PI;
    const double theta = phase + spot.longitude;
    const double phi = M_PI_2 - spot.latitude;
    center = Point(sin(phi)*cos(theta),
                   sin(phi)*sin(theta),
                   cos(phi));
    center.rotate_y(spot.star->inclination - M_PI_2);

    const double depth = sqrt(1-radius*radius);
    circle_center = {center.x*depth, center.y*depth, center.z*depth};

    const double a_x = 0;
    const double a_y = -circle_center.z/sqrt(circle_center.y*circle_center.y + circle_center.z*circle_center.z);
    const double a_z = sqrt(1-a_y*a_y);
    a = {a_x, a_y, a_z};

    const double b_x = circle_center.y*a_z - circle_center.z*a_y;
    const double b_y = circle_center.z*a_x - circle_center.x*a_z;
    const double b_z = circle_center.x*a_y - circle_center.y*a_x;
    b = {b_x, b_y, b_z};

    // Look for the smallest x value that the boundary touches to determine if the spot is entirely on the visible hemisphere
    const double theta_x_max = -2 * atan((a.x - sqrt(a.x * a.x + b.x * b.x)) / b.x);
    const double theta_x_min = -2 * atan((a.x + sqrt(a.x * a.x + b.x * b.x)) / b.x);

    const double x1 = circle_center.x + radius*(cos(theta_x_max)*a.x + sin(theta_x_max)*b.x);
    const double x2 = circle_center.x + radius*(cos(theta_x_min)*a.x + sin(theta_x_min)*b.x);
    
    is_on_edge = x1 < 0 || x2 < 0;
    visible = x1 > 0 || x2 > 0;

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


bool on_star(const double y, const double z) {
    return (y*y + z*z) <= 1;
}


bool BoundingShape::on_spot(const double y, const double z) const {
    if (!on_star(y, z)) {
        return false;
    }
    double x = sqrt(1-(y*y + z*z));
    double distance_squared = (y-circle_center.y)*(y-circle_center.y) + (z-circle_center.z)*(z-circle_center.z) + (x-circle_center.x)*(x-circle_center.x);
    return (distance_squared <= (radius*radius));
}


bounds BoundingShape::z_bounds(const double y) const {
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

    // Ignore the edge case where the brute force misses the spot entirely
    if (z_min > 1.0 || z_max > 1.0) {
        z_min = 0.0;
        z_max = 0.0;
    }

    return {z_min, z_max};
}


bounds BoundingShape::z_bounds_edge2(const double y) const {
    const double a = circle_center.z*circle_center.z + circle_center.x*circle_center.x;
    const double b = -2.0*circle_center.z*((1 + circle_center.x*circle_center.x + circle_center.y*circle_center.y + circle_center.z*circle_center.z - radius*radius)/2.0 - y*circle_center.y);
    const double c = std::pow((1 + circle_center.x*circle_center.x + circle_center.y*circle_center.y + circle_center.z*circle_center.z - radius*radius)/2.0 - y*circle_center.y, 2) *  - circle_center.x*circle_center.x *(1 - y*y) ;

    const double z_min = (-b+sqrt(b*b-4*a*c))/(2.0*a);
    const double z_max = (-b-sqrt(b*b-4*a*c))/(2.0*a);
    return {z_min, z_max};
}


bool BoundingShape::collides_with(const BoundingShape other) const {
    double distance = sqrt(
            std::pow(center.x-other.center.x, 2) +
            std::pow(center.y-other.center.y, 2) +
            std::pow(center.z-other.center.z, 2)
    );
    return distance < (max_radius+other.max_radius);
}