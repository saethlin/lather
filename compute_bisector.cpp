#include "compute_bisector.hpp"


#include <iostream>
std::vector<double> compute_bisector(const std::vector<double>& rv, const std::vector<double>& profile) {
    std::vector<double> bisector(1000);

    // Clip off the wiggles in the tails of the profile
    auto peak = std::min_element(profile.begin(), profile.end());

    auto right_end = peak;
    for (right_end = peak+1; right_end != profile.end(); right_end++) {
        if (*right_end < *(right_end-1)) {
            break;
        }
    }

    auto left_start = peak;
    for (left_start = peak-1; left_start != profile.begin(); left_start--) {
        if (*left_start < *(left_start+1)) {
            break;
        }
    }
    left_start++;

    auto left_size = peak - left_start;
    auto right_size = right_end - peak;
    if (left_size > right_size) {
        left_start = peak - right_size;
    }
    else {
        right_end = peak + left_size;
    }

    std::vector<double> left_profile(left_start, peak);
    std::vector<double> right_profile(peak, right_end);

    std::vector<double> left_rv(rv.begin() + (left_start - profile.begin()), rv.begin() + (peak - profile.begin()));
    std::vector<double> right_rv(rv.begin() + (peak - profile.begin()), rv.begin() + (right_end - profile.begin()));

    // Reverse the left side so that the profile values (x coordinate for the interpolator) are increasing
    std::reverse(left_profile.begin(), left_profile.end());
    std::reverse(left_rv.begin(), left_rv.end());

    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, left_profile.size());

    gsl_spline_init(spline, left_profile.data(), left_rv.data(), left_profile.size());

    auto interp_range = std::fabs(left_profile[0] - left_profile[left_profile.size()-1]);
    auto interp_step = interp_range / bisector.size();
    for (int i = 0; i < bisector.size(); i++ ) {
        bisector[i] = gsl_spline_eval(spline, left_profile[0] + i*interp_step, acc);
    }

    // Reverse the bisector so things are in the right order after the previous reversal
    //std::reverse(bisector.begin(), bisector.end());

    // Interpolate the right side and average
    gsl_spline* right_spline = gsl_spline_alloc(gsl_interp_cspline, right_profile.size());
    gsl_spline_init(right_spline, right_profile.data(), right_rv.data(), right_profile.size());

    // TODO Figure out why this needs to be a 1 and document it
    for (int i = 1; i < bisector.size(); i++ ) {
        bisector[i] += gsl_spline_eval(spline, right_profile[0] + i*interp_step, acc);
        bisector[i] /= 2;
    }

    gsl_spline_free(spline);
    gsl_spline_free(right_spline);
    gsl_interp_accel_free(acc);

    return bisector;
}
