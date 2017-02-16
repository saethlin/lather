#include "star.hpp"


const double solar_radius = 696000.0;
const double c = 299792458;
const double h = 6.62606896e-34;
const double k_b = 1.380e-23;
const double days_to_seconds = 86400;


double planck(double wavelength, double temperature) {
    return 2*h*c*c*1./(wavelength*wavelength*wavelength*wavelength*wavelength*(expm1((h*c)/(wavelength*k_b*temperature))));
}


Star::Star() {}


Star::Star(double radius, double period, double inclination, double temperature, double spotTempDiff,
           double limbLinear, double limbQuadratic, size_t gridSize) {
    this -> inclination = inclination * M_PI/180.0;
    this -> period = period;
    double vrot = (2*M_PI * radius*solar_radius)/(period * days_to_seconds);
    this -> equatorial_velocity = vrot * sin(this->inclination);
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

    profileQuiet = Profile(rv, ccfQuiet);
    profileActive = Profile(rv, ccfActive);

    std::ostringstream pathstream;
    pathstream << ".lathercache" << "_grid" << gridSize << "_incl" << inclination << "_period" << period;

    auto cache_path = pathstream.str();
    std::ifstream cacheread(cache_path, std::ifstream::in);

    if (cacheread.good()) {
        cacheread >> fluxQuiet;
        while (cacheread.good()) {
            cacheread >> num;
            integrated_ccf.push_back(num);
        }
        cacheread.close();
    }

    else {
        integrated_ccf = std::vector<double>(profileQuiet.size());

        for (auto y = -1.0; y <= 1; y += grid_interval) {
            std::cout << y << std::endl;
            auto& ccfShifted = quiet_profile(y);

            double limb_integral = 0.0;
            double z_bound = sqrt(1-y*y);
            for (auto z = -z_bound; z <= z_bound; z+= grid_interval) {
                double x = sqrt(1 - y*y - z*z);
                limb_integral += limb_brightness(x);
            }

            //std::cout << limb_integral << '\n';

            for (auto k = 0; k < profileQuiet.size(); k++) {
                integrated_ccf[k] += ccfShifted[k] * limb_integral;
            }

            fluxQuiet += limb_integral;
        }

        // Write out to the cache
        std::ofstream cachewrite(cache_path, std::ofstream::out);
        cachewrite << fluxQuiet << '\n' << '\n';
        for (const auto& num : integrated_ccf) {
            cachewrite << num << '\n';
        }
        cachewrite.close();
    }

    // Compute the rv that will be fitted with no spots visible.
    std::vector<double> normProfile(integrated_ccf);
    normalize(normProfile);

    fit_result = std::vector<double>(4);
    fit_result[0] = normProfile[normProfile.size()/2] - normProfile[0];
    fit_result[1] = profileQuiet.rv()[normProfile.size()/2];
    fit_result[2] = 2.71; // TODO Remove this magic number
    fit_result[3] = normProfile[0];

    fit_rv(profileQuiet.rv(), normProfile, fit_result);
    zero_rv = fit_result[1];

    //for (const auto& val : normProfile) std::cout << val << '\n';
    //exit(0);
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


double Star::limb_integral(double z_upper, double z_lower, const double y) const {
    if (z_upper < z_lower) {
        auto tmp = z_upper;
        z_upper = z_lower;
        z_lower = tmp;
    }
    bool sign = z_upper >= 0;
    double limb_integral = -1./6. * z_upper * (2 * (limbQuadratic * (3*y*y + z_upper*z_upper - 6) - 3) -3*limbLinear*(-2)) -
              1./2. * (y*y - 1.) * (limbLinear - 2*limbQuadratic) * sign*M_PI_2;

    sign = z_lower >= 0;
    limb_integral -= -1./6. * z_lower * (2 * (limbQuadratic * (3*y*y + z_lower*z_lower - 6) - 3) -3*limbLinear*(-2)) -
               1./2. * (y*y - 1.) * (limbLinear - 2*limbQuadratic) * sign*M_PI_2;

    return limb_integral;
}
