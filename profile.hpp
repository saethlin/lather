#ifndef PROFILE_HPP
#define PROFILE_HPP

#include <vector>
#include <cstddef>

class Profile {
public:
    Profile() {}
    Profile(const std::vector<double> rv, const std::vector<double> ccf, const double v_max, const double grid_interval);
    std::vector<double>& shift(const double v_shift) const;
    size_t size() const {return rv().size();}
    const std::vector<double>& rv() const {return rv_impl;}
    const std::vector<double>& ccf() const {return ccf_impl;}

private:
    double stepsize, v_max, grid_interval;
    std::vector<double> rv_impl;
    std::vector<double> ccf_impl;
    std::vector<double> derivative;
    mutable std::vector<std::vector<double> > cache;
};


#endif