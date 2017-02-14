#include "simulation.hpp"

const double pi = M_PI;


void normalize(std::vector<double>& vec) {
    auto max = *std::max_element(vec.begin(), vec.end());
    for (auto& elem : vec) {
        elem /= max;
    }
}


Simulation::Simulation(size_t gridSize, size_t spotResolution) {
    this -> gridSize = gridSize;
    this -> spotResolution = spotResolution;
}


Simulation::Simulation(const char* filename) {

    INIReader reader(filename);
    if (reader.ParseError() < 0) {
        throw std::exception();
    }

    gridSize = reader.GetInteger("simulation", "grid_resolution", 100);
    spotResolution = reader.GetInteger("simulation", "spot_resolution", 100);

    double radius = reader.GetReal("star", "radius", 1.0);
    double period = reader.GetReal("star", "period", 25.05);
    double inclination = reader.GetReal("star", "inclination", 90.0);
    double temperature = reader.GetReal("star", "Tstar", 5778.0);
    double spot_temp_diff = reader.GetReal("star", "Tdiff_spot", 663.0);
    double limbLinear = reader.GetReal("star", "limb1", 0.29);
    double limbQuadratic = reader.GetReal("star", "limb2", 0.34);

    setStar(radius, period, inclination, temperature, spot_temp_diff, limbLinear, limbQuadratic);

    for (const auto& section : reader.GetSections()) {
        if (section.substr(0, 4) == "spot") {
            double latitude = reader.GetReal(section, "latitude", 0.0);
            double longitude = reader.GetReal(section, "longitude", 180.0);
            double size = reader.GetReal(section, "size", 0.1);
            bool plage = reader.GetBoolean(section, "plage", false);

            addSpot(latitude, longitude, size, plage);
        }
    }
}


void Simulation::setStar(double radius, double period, double inclination, double temperature, double spot_temp_diff, double linear_limb, double quadratic_limb) {
    star = Star(radius, period, inclination, temperature, spot_temp_diff, linear_limb, quadratic_limb, gridSize);
}


void Simulation::addSpot(double latitude, double longitude, double fillfactor, bool plage) {
    spots.emplace_back(star, latitude, longitude, fillfactor, plage, spotResolution);
}


void Simulation::clear_spots() {
    this->spots.clear();
}



std::vector<double> Simulation::observe_rv(std::vector<double>& time, double wavelength) {
    std::vector<double> rv(time.size());

    star.intensity = planck(wavelength, star.temperature);
    for (auto& spot : spots) {
        spot.intensity = planck(wavelength, spot.temperature) / star.intensity;
    }


    std::vector<double> spot_profile(star.profileActive.size());
    std::vector<double> fit_result = star.fit_result;

    for (auto t = 0; t < time.size(); t++) {
        auto phase = fmod(time[t], star.period)/star.period * 2*pi;
        for (auto& spot : spots) {
            auto profile = spot.get_ccf(phase, wavelength);
            for (auto i = 0; i < profile.size(); i++) {
                spot_profile[i] += profile[i];
            }
        }

        // Compute the observed profile and fit the rv: the star's quiet profile minus the spot flux
        for (auto i = 0; i < spot_profile.size(); i++) {
            spot_profile[i] = star.integrated_ccf[i] - spot_profile[i];
        }
        normalize(spot_profile);
        fit_rv(star.profileQuiet.rv(), spot_profile, fit_result);
        rv[t] = fit_result[1] - star.zero_rv;

        for (auto& elem : spot_profile) elem = 0.0;
    }
    return rv;
}


std::vector<double> Simulation::observe_flux(std::vector<double>& time, double wavelength) {
    std::vector<double> flux(time.size());
    star.intensity = planck(wavelength, star.temperature);

    for (auto& spot : spots) {
        spot.intensity = planck(wavelength, spot.temperature) / star.intensity;
    }

    for (auto t = 0; t < time.size(); t++) {
        auto phase = fmod(time[t], star.period)/star.period * 2*pi;
        double spot_flux = 0.0;

        for (auto& spot : spots) {
            spot_flux += spot.get_flux(phase, wavelength);
        }

        flux[t] = (star.fluxQuiet - spot_flux) / star.fluxQuiet;
    }
    normalize(flux);
    return flux;
}


/*
void Simulation::fit(std::vector<double>& time, std::vector<double>& flux) {
    fit_sim(this, time, flux);
}
*/