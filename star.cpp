#include "star.hpp"


const double pi = M_PI;
const double solarRadius = 696000.0;

const double c = 299792458;
const double h = 6.62606896e-34;
const double k_b = 1.380e-23;
const double days_to_seconds  = 86400;


double planck(double wavelength, double temperature) {
    return 2*h*c*c*1./(wavelength*wavelength*wavelength*wavelength*wavelength*(expm1((h*c)/(wavelength*k_b*temperature))));
}


Star::Star() {}


Star::Star(double radius, double period, double inclination, double temperature, double spotTempDiff,
           double limbLinear, double limbQuadratic, size_t gridSize) {
    this -> period = period;
    this -> vrot = (2*pi * radius*solarRadius)/(period * days_to_seconds);
    this -> inclination = inclination * pi/180.0;
    this -> temperature = temperature;
    this -> spotTempDiff = spotTempDiff;
    this -> limbLinear = limbLinear;
    this -> limbQuadratic = limbQuadratic;
    this -> gridSize = gridSize;
    this -> intensity = planck(5293.4115e-10, temperature);

    std::vector<double> grid_steps(gridSize+1);
    for (auto i = 0; i <= gridSize; i++) {
        grid_steps[i] = -1.0 + i * 2.0 / gridSize;
    }

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
    pathstream << ".lathercache";
    pathstream << "_grid" << gridSize;
    pathstream << "_incl" << inclination;
    pathstream << "_period" << period;

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

        for (const auto& y : grid_steps) {
            auto v_shift = y * vrot * sin(this->inclination);
            auto& ccfShifted = profileQuiet.shift(v_shift);
            auto limbSum = 0.0;

            for (const auto& z : grid_steps) {
                auto rSquared = y*y + z*z;

                if (rSquared <= 1) {
                    auto r_cos = sqrt(1 - rSquared);
                    limbSum += 1 - limbLinear * (1. - r_cos) - limbQuadratic * (1 - r_cos) * (1 - r_cos);
                }
            }

            for (auto k = 0; k < profileQuiet.size(); k++) {
                integrated_ccf[k] += ccfShifted[k] * limbSum;
            }

            fluxQuiet += limbSum;
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
}
