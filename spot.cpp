#include "spot.hpp"


Spot::Spot(Star* star, double latitude, double longitude, double fillfactor, bool plage) {
    this->latitude = latitude * M_PI/180.0;
    this->longitude = longitude * M_PI/180.0;
    size = sqrt(2*fillfactor);
    this->star = star;
    this->plage = plage;
    this->temperature = star->temperature - star->spotTempDiff;
    auto h = 1-sqrt(1-size*size);
    radius_squared = 2*h - h*h + h*h;
}


spot_bounds Spot::get_bounds_at(double phase) {
    double theta = phase + longitude;
    double phi = M_PI_2 - latitude;
    center = Point(sin(phi)*cos(theta),
                   sin(phi)*sin(theta),
                   cos(phi));
    center.rotate_y(star->inclination - M_PI_2);

    bool visible = center.x > -sqrt(2*size);

    double depth = sqrt(1-size*size);
    Point circle_center(center.x*depth, center.y*depth, center.z*depth);

    auto h = 1-sqrt(1-size*size);
    auto r = sqrt(2*h - h*h);

    double a_x = 0;
    double a_y = -circle_center.z/sqrt(circle_center.y*circle_center.y + circle_center.z*circle_center.z);
    double a_z = sqrt(1-a_y*a_y);

    double b_x = circle_center.y*a_z - circle_center.z*a_y;
    double b_y = circle_center.z*a_x - circle_center.x*a_z;
    double b_z = circle_center.x*a_y - circle_center.y*a_x;

    double theta_y_max = -2 * atan(a_y/b_y - (sqrt(a_y*a_y + b_y*b_y)/b_y));
    double theta_y_min = -2 * atan(a_y/b_y + (sqrt(a_y*a_y + b_y*b_y)/b_y));

    double y_max = circle_center.y + r*cos(theta_y_max)*a_y + r*sin(theta_y_max)*b_y;
    double y_min = circle_center.y + r*cos(theta_y_min)*a_y + r*sin(theta_y_min)*b_y;

    double theta_z_max = -2 * atan(a_z/b_z - (sqrt(a_z*a_z + b_z*b_z)/b_z));
    double theta_z_min = -2 * atan(a_z/b_z + (sqrt(a_z*a_z + b_z*b_z)/b_z));

    double z_max = circle_center.z + r*cos(theta_z_max)*a_z + r*sin(theta_z_max)*b_z;
    double z_min = circle_center.z + r*cos(theta_z_min)*a_z + r*sin(theta_z_min)*b_z;

    y_min = floor(y_min/star->grid_interval)*star->grid_interval;
    y_max = ceil(y_max/star->grid_interval)*star->grid_interval;
    z_min = floor(z_min/star->grid_interval)*star->grid_interval;
    z_max = ceil(z_max/star->grid_interval)*star->grid_interval;

    return {visible, y_min, y_max, z_min, z_max};
}


std::vector<double> Spot::get_ccf(double phase, double wavelength) {
    std::vector<double> profile(star->profileActive.size());
    auto bounds = get_bounds_at(phase);
    if (not bounds.visible) {
        return profile;
    }

    for (double y = bounds.miny; y < bounds.maxy; y += star->grid_interval) {

        auto& ccfShifted = star->quiet_profile(y);
        auto& ccfActiveShifted = star->active_profile(y);
        double limb_integral = 0.0;

        for (double z = bounds.minz; z < bounds.maxz; z += star->grid_interval) {
            auto x = sqrt(1 - z*z - y*y);
            if (is_on_spot(x, y, z)) {
                limb_integral += star->limb_brightness(x);
            }
        }

        for (auto i = 0; i < ccfShifted.size(); i++) {
            profile[i] += (ccfShifted[i] - intensity * ccfActiveShifted[i]) * limb_integral;
        }
    }
    return profile;
}


double Spot::get_flux(double phase, double wavelength) {
    auto bounds = get_bounds_at(phase);
    if (not bounds.visible) {
        return 0.0;
    }

    double limb_integral = 0.0;
    for (auto y = bounds.miny; y < bounds.maxy; y += star->grid_interval) {
        for (auto z = bounds.minz; z < bounds.maxz; z += star->grid_interval) {
            auto x = sqrt(1 - y*y - z*z);
            if (is_on_spot(x, y, z)) {
                limb_integral += star->limb_brightness(x);
            }
        }
    }
    auto intensity = planck(wavelength, temperature) / star->intensity;
    return (1-intensity)*limb_integral;
}


bool Spot::is_on_spot(double x, double y, double z) {
    if ((y*y + z*z) >= 1.0) { // Check if on the star disk
        return false;
    }
    else {
        double distance_squared = (x-center.x)*(x-center.x) + (y-center.y)*(y-center.y) + (z-center.z)*(z-center.z);
        return distance_squared < radius_squared;
    }
}
