#ifndef PROFILE_HPP
#define PROFILE_HPP
#include <map>
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
        std::map<double, std::vector<double>* > cache;
};
#endif