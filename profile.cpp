#include <vector>
#include <map>
#include <math.h>
#include "profile.hpp"


Profile::Profile() {}


Profile::Profile(std::vector<double> rv, std::vector<double> ccf) {
    this -> rv = rv;
    this -> ccf = ccf;
    this -> size = rv.size();

    cache = std::map<double, std::vector<double>* > ();

    derivative = std::vector<double>(size);
    for (unsigned int i = 0; i < size-1; i++) {
        derivative[i] = (ccf[i+1] - ccf[i]) / (rv[i+1] - rv[i]);
    }
    stepsize = rv[1]-rv[0];
}


std::vector<double> Profile::shift(double v_shift) {

    if (!cache.count(v_shift)) {

        cache[v_shift] = new std::vector<double> (size);
        std::vector<double> &ccf_shifted = *cache[v_shift];

        unsigned int i;
        int quotient = round(v_shift / stepsize);
        double remainder = v_shift - ((double)quotient)*stepsize;


        if (v_shift >= 0) {
            for (i = 0; i < (unsigned int)(quotient+1); i++) {
                ccf_shifted[i] = ccf[0];
            }
            for (i = quotient+1; i < size; i++) {
                ccf_shifted[i] = ccf[i-quotient] - remainder * derivative[i-quotient];
            }
        }
        else {
            for (i = 1; i < size+quotient; i++) {
                ccf_shifted[i] = ccf[i - quotient] - remainder * derivative[i - quotient - 1];
            }
            for (i = size+quotient; i < size; i++) {
                ccf_shifted[i] = ccf[size-1];
            }
            ccf_shifted[0] = ccf[0];
        }
    }
    return *cache[v_shift];
}
