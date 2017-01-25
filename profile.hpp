#ifndef PROFILE_HPP
#define PROFILE_HPP
#include <cmath>
#include <unordered_map>
#include <memory>
#include <vector>

class Profile {
    public:
        Profile();

        Profile(std::vector<double> rv, std::vector<double> ccf);

        unsigned int size;
        std::vector<double> rv;
        std::vector<double> ccf;

        std::vector<double> shift(double v_shift);

    private:
        double stepsize;
        std::vector<double> derivative;
        std::unordered_map<double, std::shared_ptr<std::vector<double> > > cache;
};
#endif