#include "spot.hpp"


Spot::Spot(Star* star, double latitude, double longitude, double fillfactor, bool plage) {
    this->latitude = latitude * M_PI/180.0;
    this->longitude = longitude * M_PI/180.0;
    size = sqrt(2*fillfactor);
    this->star = star;
    this->plage = plage;
    this->temperature = star->temperature - star->spotTempDiff;
}


double Spot::get_flux(double phase) {
    double limb_integral = 0.0;
    auto bounds = BoundingShape(*this, phase);
    if (not bounds.is_visible()) {
        return limb_integral;
    }

    auto y_bounds = bounds.get_y_bounds();
    for (auto y = y_bounds.lower; y < y_bounds.upper; y += star->grid_interval) {
        auto z_bounds = bounds.get_z_bounds(y);
        limb_integral += star->get_limb_integral(y, z_bounds.lower, z_bounds.upper);
    }
    return (1-intensity)*limb_integral;
}


std::vector<double> Spot::get_ccf(double time) {
    std::vector<double> profile(star->profileActive.size());
    auto bounds = BoundingShape(*this, time);
    if (not bounds.is_visible()) {
        return profile;
    }

    auto y_bounds = bounds.get_y_bounds();
    for (double y = y_bounds.lower; y < y_bounds.upper; y += star->grid_interval) {

        auto& ccfShifted = star->quiet_profile(y);
        auto& ccfActiveShifted = star->active_profile(y);

        auto z_bounds = bounds.get_z_bounds(y);
        double limb_integral = star->get_limb_integral(z_bounds.lower, z_bounds.upper, y);

        for (auto i = 0; i < ccfShifted.size(); i++) {
            profile[i] += (ccfShifted[i] - intensity * ccfActiveShifted[i]) * limb_integral;
        }
    }

    return profile;
}
