#include "star.hpp"


const double solar_radius = 696000.0;
const double days_to_seconds = 86400;


Star::Star(double radius, double period, double inclination, double temperature, double spotTempDiff,
           double limbLinear, double limbQuadratic, size_t gridSize) {
    this -> inclination = inclination * M_PI/180.0;
    this -> period = period;
    double edge_velocity = (2*M_PI * radius*solar_radius)/(period * days_to_seconds);
    this -> equatorial_velocity = edge_velocity * sin(this->inclination);
    this -> temperature = temperature;
    this -> spotTempDiff = spotTempDiff;
    this -> limbLinear = limbLinear;
    this -> limbQuadratic = limbQuadratic;
    this -> grid_interval = 2.0/gridSize;

    // Setup for profiles
    std::vector<double> rv;
    std::vector<double> ccfQuiet;
    std::vector<double> ccfActive;

    // SOAP makes a decision between low or high res ccf based on rotational velocity
    std::ifstream ifs("/home/ben/lather/resources/solarccfhires.txt");
    std::string line;
    double num;

    // Skip first two header lines
    getline(ifs, line);
    getline(ifs, line);

    while (ifs.good()) {
        ifs >> num;
        rv.push_back(num);

        ifs >> num;
        ccfQuiet.push_back(num);

        ifs >> num;
        ccfActive.push_back(num);
    }
    ifs.close();

    profileQuiet = Profile(rv, ccfQuiet, equatorial_velocity, grid_interval);
    profileActive = Profile(rv, ccfActive, equatorial_velocity, grid_interval);

    integrated_ccf = std::vector<double>(profileQuiet.size());

    for (double y = -1.0; y <= 1; y += grid_interval) {

        auto& ccfShifted = quiet_profile(y);
        double z_bound = sqrt(1 - y*y);
        double limb_integral = get_limb_integral(z_bound, -z_bound, y);

        for (auto k = 0; k < profileQuiet.size(); k++) {
            integrated_ccf[k] += ccfShifted[k] * limb_integral;
        }

        fluxQuiet += limb_integral;
    }

    // Compute the rv that will be fitted with no spots visible.
    std::vector<double> normProfile(integrated_ccf);
    normalize(normProfile);

    fit_result = std::vector<double>(4);
    fit_result[0] = normProfile[normProfile.size()/2] - normProfile[0];
    fit_result[1] = profileQuiet.rv()[normProfile.size()/2];
    fit_result[2] = 2.71; // TODO Remove this magic number
    fit_result[3] = normProfile[0];

    fit_result = fit_rv(profileQuiet.rv(), normProfile, fit_result);
    zero_rv = fit_result[1];
}


double Star::limb_brightness(const double r_cos) const {
    return 1 - limbLinear * (1 - r_cos) - limbQuadratic * (1 - r_cos) * (1 - r_cos);
}


std::vector<double>& Star::quiet_profile(const double y) {
    return profileQuiet.shift(y * equatorial_velocity);
}


std::vector<double>& Star::active_profile(const double y) {
    return profileActive.shift(y * equatorial_velocity);
}


double Star::get_limb_integral(const double z_upper, const double z_lower, const double y) const {
    if (z_lower == z_upper) return 0;

    double x_upper = sqrt(1 - std::min(z_upper*z_upper + y*y, 1.0));
    double x_lower = sqrt(1 - std::min(z_lower*z_lower + y*y, 1.0));

    double limb_integral = 1./6. * (z_upper * (3*limbLinear*(x_upper-2) + 2*(limbQuadratic*(3*x_upper + 3*y*y + z_upper*z_upper - 6) + 3 )) -
            3 * (y*y - 1)*(limbLinear + 2*limbQuadratic)*atan(z_upper/x_upper));

    limb_integral -= 1./6. * (z_lower * (3*limbLinear*(x_lower-2) + 2*(limbQuadratic*(3*x_lower + 3*y*y + z_lower*z_lower - 6) + 3 )) -
                                    3 * (y*y - 1)*(limbLinear + 2*limbQuadratic)*atan(z_lower/x_lower));

    return limb_integral;
}
