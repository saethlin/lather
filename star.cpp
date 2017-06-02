#include "star.hpp"


const double solar_radius = 696000.0;
const double days_to_seconds = 86400;


Star::Star(const double radius, const double period, const double inclination, const double temperature,
           const double spot_temp_diff, const double limb_linear, const double limb_quadratic, size_t gridSize) {
    this->inclination = inclination * M_PI / 180.0;
    this->period = period;
    double edge_velocity = (2 * M_PI * radius * solar_radius) / (period * days_to_seconds);
    this->equatorial_velocity = edge_velocity * sin(this->inclination);
    this->temperature = temperature;
    this->spot_temp_diff = spot_temp_diff;
    this->limb_linear = limb_linear;
    this->limb_quadratic = limb_quadratic;
    this->grid_interval = 2.0 / gridSize;

    // Seed the random number generator
    generator.seed(0);

    // Setup for profiles
    std::vector<double> rv;
    std::vector<double> ccf_quiet;
    std::vector<double> ccf_active;
    {
        std::ifstream ifs("/home/ben/lather/resources/solarccfhires.txt");
        std::string line;

        // Skip first two header lines
        getline(ifs, line);
        getline(ifs, line);

        double rv_val, quiet_val, active_val;
        while (ifs >> rv_val >> quiet_val >> active_val) {
            rv.push_back(rv_val);
            ccf_quiet.push_back(quiet_val);
            ccf_active.push_back(active_val);
        }
    }

    profile_quiet = Profile(rv, ccf_quiet, equatorial_velocity, grid_interval);
    profile_active = Profile(rv, ccf_active, equatorial_velocity, grid_interval);

    integrated_ccf = std::vector<double>(profile_quiet.size());

    for (double y = -1.0; y <= 1.0; y += grid_interval) {

        auto& ccfShifted = quiet_profile(y);
        double z_bound = sqrt(1 - y*y);
        double limb_integral = get_limb_integral(-z_bound, z_bound, y);

        for (auto k = 0; k < profile_quiet.size(); k++) {
            integrated_ccf[k] += ccfShifted[k] * limb_integral;
        }

        flux_quiet += limb_integral;
    }

    // Compute the rv that will be fitted with no spots visible.
    std::vector<double> normalized_profile(integrated_ccf);
    normalize(normalized_profile);

    fit_result = std::vector<double>(4);
    fit_result[0] = normalized_profile[normalized_profile.size()/2] - normalized_profile[0];
    fit_result[1] = profile_quiet.rv()[normalized_profile.size()/2];
    fit_result[2] = 2.71; // TODO Remove this magic number
    fit_result[3] = normalized_profile[0];

    fit_result = fit_rv(profile_quiet.rv(), normalized_profile, fit_result);
    zero_rv = fit_result[1];

    // Draw the star
    image.reserve(gridSize*gridSize*4);
    for (double y = -1.0; y <= 1.0; y += grid_interval) {
        for (double z = -1.0; z <= 1.0; z += grid_interval) {
            float intensity = 0.0;
            if (y*y + z*z <= 1.0) {
                double x = 1-(y*y+z*z);
                x = std::max(x, 0.0);
                intensity = (float)limb_brightness(x);
            }
            image.push_back(intensity*255);
            image.push_back(intensity*157);
            image.push_back(intensity*63);
            image.push_back(0.0);
        }
    }
}


double Star::limb_brightness(const double r_cos) const {
    return 1 - limb_linear * (1 - r_cos) - limb_quadratic * (1 - r_cos) * (1 - r_cos);
}


std::vector<double>& Star::quiet_profile(const double y) const {
    return profile_quiet.shift(y * equatorial_velocity);
}


std::vector<double>& Star::active_profile(const double y) const {
    return profile_active.shift(y * equatorial_velocity);
}


double Star::get_limb_integral(const double z_lower, const double z_upper, const double y) const {
    if (z_lower == z_upper) return 0;

    double x_upper = sqrt(1 - std::min(z_upper*z_upper + y*y, 1.0));
    double x_lower = sqrt(1 - std::min(z_lower*z_lower + y*y, 1.0));

    double limb_integral = 1./6. * (z_upper * (3*limb_linear*(x_upper-2) + 2*(limb_quadratic*(3*x_upper + 3*y*y + z_upper*z_upper - 6) + 3 )) -
            3 * (y*y - 1)*(limb_linear + 2*limb_quadratic)*atan(z_upper/x_upper));

    limb_integral -= 1./6. * (z_lower * (3*limb_linear*(x_lower-2) + 2*(limb_quadratic*(3*x_lower + 3*y*y + z_lower*z_lower - 6) + 3 )) -
                                    3 * (y*y - 1)*(limb_linear + 2*limb_quadratic)*atan(z_lower/x_lower));

    return limb_integral;
}
