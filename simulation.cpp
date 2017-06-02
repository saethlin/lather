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

    dynamic_fill_factor = reader.GetReal("star", "fillfactor", 0.0);

    for (const auto& section : reader.GetSections()) {
        if (section.substr(0, 4) == "spot") {
            const double latitude = reader.GetReal(section, "latitude", 0.0);
            const double longitude = reader.GetReal(section, "longitude", 180.0);
            const double size = reader.GetReal(section, "size", 0.1);
            const bool plage = reader.GetBoolean(section, "plage", false);

            add_spot(latitude, longitude, size, plage, false);
        }
    }
}


void Simulation::set_star(double radius, double period, double inclination, double temperature,
                          double spot_temp_diff, double linear_limb, double quadratic_limb, size_t grid_size) {
    star = Star(radius, period, inclination, temperature, spot_temp_diff, linear_limb, quadratic_limb, grid_size);
}


void Simulation::add_spot(double latitude, double longitude, double fillfactor, bool plage, bool mortal) {
    spots.emplace_back(&star, latitude, longitude, fillfactor, plage, mortal);
}


void Simulation::clear_spots() {
    spots.clear();
}


void Simulation::check_fill_factor(double time) {
    double current_fill_factor = 0.0;
    for (const auto& spot : spots) {
        if (spot.alive(time)) {
            current_fill_factor += (spot.radius*spot.radius)/2.0;
        }
    }

    std::lognormal_distribution<> fill_dist(0.5, 4.0);
    std::uniform_real_distribution<> lat_dist(-30.0, 30.0);
    std::uniform_real_distribution<> long_dist(0.0, 360.0);

    while (current_fill_factor < dynamic_fill_factor) {
        double new_fill_factor = fill_dist(star.generator)*9.4e-6;
        while (new_fill_factor > 0.001) {
            new_fill_factor = fill_dist(star.generator)*9.4e-6;
        }

        Spot new_spot(&star, lat_dist(star.generator), long_dist(star.generator), new_fill_factor, false, true);
        new_spot.time_disappear += time;
        new_spot.time_appear += time;

        bool collides = false;
        for (const auto& spot: spots) {
            if (spot.alive(new_spot.time_disappear) || spot.alive(new_spot.time_appear)) {
                if (new_spot.collides_with(spot)) {
                    collides = true;
                    break;
                }
            }
        }
        if (! collides) {
            current_fill_factor += (new_spot.radius*new_spot.radius)/2.0;
            spots.push_back(new_spot);
        }
    }
}


std::vector<double> Simulation::observe_rv(const std::vector<double>& time, const double wavelength_min, const double wavelength_max) {
    std::vector<double> rv(time.size());

    star.intensity = planck_integral(star.temperature, wavelength_min, wavelength_max);

    auto fit_guess = star.fit_result;

    // Ensure no collisions at every time of observation
    for (const auto& t : time) check_fill_factor(t);

    double intensity = planck_integral((star.temperature-star.spot_temp_diff), wavelength_min, wavelength_max) / star.intensity;
    for (auto &spot : spots) spot.intensity = intensity;

    for (auto t = 0; t < time.size(); t++) {
        check_fill_factor(time[t]);
        std::vector<double> spot_profile(star.profile_active.size());
        for (const auto &spot : spots) {
            if (spot.alive(time[t])) {
                auto profile = spot.get_ccf(time[t]);
                for (auto i = 0; i < profile.size(); i++) {
                    spot_profile[i] += profile[i];
                }
            }
        }

        // Compute the observed profile and fit the rv: the star's quiet profile minus the spot flux
        for (auto i = 0; i < spot_profile.size(); i++) {
            spot_profile[i] = star.integrated_ccf[i] - spot_profile[i];
        }
        normalize(spot_profile);
        auto fit_result = fit_rv(star.profile_quiet.rv(), spot_profile, fit_guess);
        rv[t] = fit_result[1] - star.zero_rv;

        // TODO: Compute profile bisector
        compute_bisector(star.profile_quiet.rv(), spot_profile);
        exit(0);
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
        check_fill_factor(time[t]);
        double spot_flux = 0.0;
        for (const auto& spot : spots) {
            spot_flux += spot.get_flux(time[t]);
        }

        flux[t] = (star.flux_quiet - spot_flux) / star.flux_quiet;
    }
    normalize(flux);
    return flux;
}


void Simulation::draw(const double time, const int i) const {
    std::vector<double> image(star.image.begin(), star.image.end());
    for (const auto& spot : spots) {
        if (spot.alive(time)) {

            const auto bounds = BoundingShape(spot, time);
            auto y_bounds = bounds.y_bounds();
            y_bounds = {round(y_bounds.upper / star.grid_interval) * star.grid_interval,
                        round(y_bounds.lower / star.grid_interval) * star.grid_interval};
            for (auto y = y_bounds.lower; y < y_bounds.upper; y += star.grid_interval) {
                auto z_bounds = bounds.z_bounds(y);
                z_bounds = {round(z_bounds.upper / star.grid_interval) * star.grid_interval,
                            round(z_bounds.lower / star.grid_interval) * star.grid_interval};
                for (auto z = z_bounds.lower; z < z_bounds.upper; z += star.grid_interval) {
                    unsigned int y_index = round((y + 1.0) / 2 * 1000);
                    unsigned int z_index = round((z + 1.0) / 2 * 1000);
                    double x = 1 - (y * y + z * z);
                    x = std::max(x, 0.0);
                    image[z_index * 1000 + y_index] = star.limb_brightness(x) * spot.intensity;
                }
            }
        }
    }
    normalize(image);
    for (auto& val : image) val = 1-val;
    Magick::InitializeMagick(".");
    Magick::Image png(1000, 1000, "K", Magick::DoublePixel, image.data());
    std::stringstream ss;
    ss << "frame" << std::setw(5) << std::setfill('0') << i << ".png";
    png.write(ss.str());
}


std::vector<uint8_t> Simulation::draw_rgba(const double time) {
    check_fill_factor(time);
    star.intensity = planck_integral(star.temperature, 4000e-10, 7000e-10);
    for (auto &spot : spots) {
        spot.intensity = planck_integral(spot.temperature, 4000e-10, 7000e-10) / star.intensity;
    }

    std::vector<float> image = star.image;
    for (const auto& spot : spots) {
        if (spot.alive(time)) {

            const auto bounds = BoundingShape(spot, time);
            auto y_bounds = bounds.y_bounds();
            y_bounds = {round(y_bounds.upper / star.grid_interval) * star.grid_interval,
                        round(y_bounds.lower / star.grid_interval) * star.grid_interval};
            for (auto y = y_bounds.lower; y < y_bounds.upper; y += star.grid_interval) {
                auto z_bounds = bounds.z_bounds(y);
                z_bounds = {round(z_bounds.upper / star.grid_interval) * star.grid_interval,
                            round(z_bounds.lower / star.grid_interval) * star.grid_interval};
                for (auto z = z_bounds.lower; z < z_bounds.upper; z += star.grid_interval) {
                    unsigned int y_index = round((y + 1.0) / 2 * 1000);
                    unsigned int z_index = round((z + 1.0) / 2 * 1000);
                    double x = 1 - (y * y + z * z);
                    x = std::max(x, 0.0);
                    float intensity = (float)(star.limb_brightness(x)*spot.intensity);
                    int index = z_index*1000 + y_index;
                    image[4*index] = intensity*255;
                    image[4*index+1] = intensity*131;
                    image[4*index+2] = 0.0;
                    image[4*index+3] = 0.0;
                }
            }
        }
    }

    // Normalize
    float max_val = *std::max_element(image.begin(), image.end());
    for (int i = 0; i < image.size()/4; i++) {
        image[4*i] *= 255/max_val;
        image[4*i+1] *= 255/max_val;
        image[4*i+2] *= 255/max_val;
        image[4*i+3] = 255;
    }
    return std::vector<uint8_t>(image.begin(), image.end());
}