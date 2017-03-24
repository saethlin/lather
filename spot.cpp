#include "spot.hpp"


Spot::Spot(Star* star, const double latitude, const double longitude, const double fillfactor, const bool plage, const bool mortal) {
    this->latitude = latitude * M_PI/180.0;
    this->longitude = longitude * M_PI/180.0;
    radius = sqrt(2*fillfactor);
    this->star = star;
    this->plage = plage;
    this->mortal = mortal;
    if (plage) {
        this->temperature = star->temperature + star->spot_temp_diff;
    }
    else {
        this->temperature = star->temperature - star->spot_temp_diff;
    }

    if (mortal) {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<> d(0, 5.0);

        time_appear = 0.0;
        time_disappear = 15.0 + d(gen);
    }
}


double Spot::get_flux(const double time) const {
    double limb_integral = 0.0;
    const auto bounds = BoundingShape(*this, time);
    const auto y_bounds = bounds.y_bounds();
    for (auto y = y_bounds.lower; y < y_bounds.upper; y += star->grid_interval) {
        const auto z_bounds = bounds.z_bounds(y);
        limb_integral += star->get_limb_integral(z_bounds.lower, z_bounds.upper, y);
    }
    return (1-intensity)*limb_integral;
}


std::vector<double> Spot::get_ccf(const double time) const {
    std::vector<double> profile(star->profile_active.size());
    const auto bounds = BoundingShape(*this, time);
    const auto y_bounds = bounds.y_bounds();
    for (double y = y_bounds.lower; y < y_bounds.upper; y += star->grid_interval) {

        const auto& ccf_quiet_shifted = star->quiet_profile(y);
        const auto& ccf_active_shifted = star->active_profile(y);

        const auto z_bounds = bounds.z_bounds(y);
        const double limb_integral = star->get_limb_integral(z_bounds.lower, z_bounds.upper, y);

        for (auto i = 0; i < ccf_quiet_shifted.size(); i++) {
            profile[i] += (ccf_quiet_shifted[i] - intensity * ccf_active_shifted[i]) * limb_integral;
        }
    }

    return profile;
}


bool Spot::alive(const double time) const {
    if (!mortal) return true;
    return time >= time_appear && time <= time_disappear;
}


bool Spot::collides_with(const Spot& other) const {
    BoundingShape bounds(*this, 0);
    BoundingShape other_bounds(other, 0);
    return bounds.collides_with(other_bounds);
}
