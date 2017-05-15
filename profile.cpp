#include "profile.hpp"
#include <iostream>


Profile::Profile(const std::vector<double> rv, const std::vector<double> ccf, const double v_max, const double grid_interval) {
    this->v_max = v_max;
    this->grid_interval = grid_interval;
    rv_impl = rv;
    ccf_impl = ccf;
    stepsize = rv[1]-rv[0];

    cache.resize((size_t)(2.0/grid_interval));

    derivative = std::vector<double>(size());
    for (auto i = 0; i < size()-1; i++) {
        derivative[i] = (ccf[i+1] - ccf[i]) / (rv[i+1] - rv[i]);
    }
}


std::vector<double>& Profile::shift(const double v_shift) const {

    const int index = (int)(v_shift/v_max/grid_interval + 1.0/grid_interval);
    auto& ccf_shifted = cache[index];

    if (ccf_shifted.size() == 0) {

        ccf_shifted = std::vector<double>(size());

        const int quotient = (int)round(v_shift / stepsize);
        const double remainder = v_shift - ((double)quotient)*stepsize;

        if (v_shift >= 0) {
            for (auto i = 0; i < quotient+1; i++) {
                ccf_shifted[i] = ccf_impl[0];
            }
            for (auto i = quotient+1; i < size(); i++) {
                ccf_shifted[i] = ccf_impl[i-quotient] - remainder * derivative[i-quotient];
            }
        }
        else {
            ccf_shifted[0] = ccf_impl[0];
            for (auto i = 1; i < size()+quotient; i++) {
                ccf_shifted[i] = ccf_impl[i-quotient] - remainder * derivative[i-quotient-1];
            }
            for (auto i = size()+quotient; i < size(); i++) {
                ccf_shifted[i] = ccf_impl[size()-1];
            }
        }
    }

    return ccf_shifted;
}
