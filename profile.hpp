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
    std::vector<double>& shift(double v_shift);
    size_t size() const {return rv().size();}
    const std::vector<double>& rv() const {return rv_impl;}
    const std::vector<double>& ccf() const {return ccf_impl;}

private:
    double stepsize;
    std::vector<double> rv_impl;
    std::vector<double> ccf_impl;
    std::vector<double> derivative;
    std::unordered_map<double, std::shared_ptr<std::vector<double> > > cache;
};


#endif