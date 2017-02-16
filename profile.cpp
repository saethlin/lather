#include "profile.hpp"


Profile::Profile() {}


Profile::Profile(std::vector<double> rv, std::vector<double> ccf) {
    rv_impl = rv;
    ccf_impl = ccf;

    cache = std::unordered_map<double, std::shared_ptr<std::vector<double> > >();

    derivative = std::vector<double>(size());
    for (auto i = 0; i < size()-1; i++) {
        derivative[i] = (ccf[i+1] - ccf[i]) / (rv[i+1] - rv[i]);
    }
    stepsize = rv[1]-rv[0];
}


std::vector<double>& Profile::shift(double v_shift) {

    auto& entry = cache[v_shift];

    if (!entry) {
        cache[v_shift] = std::make_shared<std::vector<double> >(size());
        auto& ccf_shifted = *entry;

        int quotient = round(v_shift / stepsize);
        double remainder = v_shift - ((double)quotient)*stepsize;

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
    return *entry.get();
}
