#include "simulation.hpp"


void normalize(std::vector<double>& vec) {
    const auto max = *std::max_element(vec.begin(), vec.end());
    for (auto& elem : vec) {
        elem /= max;
    }
}


Simulation::Simulation(const char* filename) {

    INIReader reader(filename);
    if (reader.ParseError() < 0) {
        throw std::exception();
    }

    const size_t grid_size = (size_t)reader.GetInteger("star", "grid_resolution", 100);
    const double radius = reader.GetReal("star", "radius", 1.0);
    const double period = reader.GetReal("star", "period", 25.05);
    const double inclination = reader.GetReal("star", "inclination", 90.0);
    const double temperature = reader.GetReal("star", "Tstar", 5778.0);
    const double spot_temp_diff = reader.GetReal("star", "Tdiff_spot", 663.0);
    const double limbLinear = reader.GetReal("star", "limb1", 0.29);
    const double limbQuadratic = reader.GetReal("star", "limb2", 0.34);

    set_star(radius, period, inclination, temperature, spot_temp_diff, limbLinear, limbQuadratic, grid_size);

    for (const auto& section : reader.GetSections()) {
        if (section.substr(0, 4) == "spot") {
            const double latitude = reader.GetReal(section, "latitude", 0.0);
            const double longitude = reader.GetReal(section, "longitude", 180.0);
            const double size = reader.GetReal(section, "size", 0.1);
            const bool plage = reader.GetBoolean(section, "plage", false);

            add_spot(latitude, longitude, size, plage);
        }
    }
}


void Simulation::set_star(double radius, double period, double inclination, double temperature,
                          double spot_temp_diff, double linear_limb, double quadratic_limb, size_t grid_size) {
    star = Star(radius, period, inclination, temperature, spot_temp_diff, linear_limb, quadratic_limb, grid_size);
}


void Simulation::add_spot(double latitude, double longitude, double fillfactor, bool plage) {
    spots.emplace_back(&star, latitude, longitude, fillfactor, plage);
}


void Simulation::clear_spots() {
    spots.clear();
}

void Simulation::check_fill_factor(double time) {

}


std::vector<double> Simulation::observe_rv(const std::vector<double>& time, const double wavelength_min, const double wavelength_max) {
    std::vector<double> rv(time.size());

    star.intensity = planck_integral(star.temperature, wavelength_min, wavelength_max);
    for (auto &spot : spots) {
        spot.intensity = planck_integral(spot.temperature, wavelength_min, wavelength_max) / star.intensity;
    }

    std::vector<double> spot_profile(star.profile_active.size());
    auto fit_guess = star.fit_result;

    for (auto t = 0; t < time.size(); t++) {
        for (const auto &spot : spots) {
            auto profile = spot.get_ccf(time[t]);
            for (auto i = 0; i < profile.size(); i++) {
                spot_profile[i] += profile[i];
            }
        }
        // Compute the observed profile and fit the rv: the star's quiet profile minus the spot flux
        for (auto i = 0; i < spot_profile.size(); i++) {
            spot_profile[i] = star.integrated_ccf[i] - spot_profile[i];
        }
        normalize(spot_profile);
        auto fit_result = fit_rv(star.profile_quiet.rv(), spot_profile, fit_guess);
        rv[t] = fit_result[1] - star.zero_rv;

        for (auto &elem : spot_profile) elem = 0.0;
    }
    for (auto& val: rv) val *= 1000.0;  // Convert to m/s
    return rv;
}


std::vector<double> Simulation::observe_flux(const std::vector<double>& time, const double wavelength_min, const double wavelength_max) {
    std::vector<double> flux(time.size());

    star.intensity = planck_integral(star.temperature, wavelength_min, wavelength_max);
    for (auto &spot : spots) {
        spot.intensity = planck_integral(spot.temperature, wavelength_min, wavelength_max) / star.intensity;
    }

    for (auto t = 0; t < time.size(); t++) {
        double spot_flux = 0.0;
        for (const auto& spot : spots) {
            spot_flux += spot.get_flux(time[t]);
        }

        flux[t] = (star.flux_quiet - spot_flux) / star.flux_quiet;
    }
    normalize(flux);
    return flux;
}


void Simulation::draw(const double time) const {
}