#include "compute_bisector.hpp"
#include <algorithm>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>


std::vector<double> compute_bisector(const std::vector<double>& rv, const std::vector<double>& profile) {
    std::vector<double> bisector(10000);

    // Clip off the wiggles in the tails of the profile
    auto peak = std::min_element(profile.begin(), profile.end());

    auto right_end = peak;
    while (*(right_end+1) >= *right_end) {
        right_end += 1;
    }

    auto left_start = peak;
    while (*(left_start-1) >= *left_start) {
        left_start -= 1;
    }

    auto left_size = peak - left_start;
    auto right_size = right_end - peak;
    if (left_size > right_size) {
        left_start = peak - right_size;
    }
    else {
        right_end = peak + left_size;
    }

    // This +1 makes the left side include the peak
    std::vector<double> left_profile(left_start, peak+1);
    std::vector<double> right_profile(peak, right_end+1);

    std::vector<double> left_rv(rv.begin() + (left_start - profile.begin()), rv.begin() + (peak - profile.begin())+1);
    std::vector<double> right_rv(rv.begin() + (peak - profile.begin()), rv.begin() + (right_end - profile.begin())+1);

    std::vector<double> profile_interp_values(bisector.size());
    auto start = *peak;
    auto stop = std::min(*left_start, *right_end);
    for (int i = 0; i < bisector.size(); i++) {
        profile_interp_values[i] = start + i * ((stop-start)/bisector.size());
    }

    // Reverse the left side so that the profile values (x coordinate for the interpolator) are increasing
    std::reverse(left_profile.begin(), left_profile.end());
    std::reverse(left_rv.begin(), left_rv.end());

    gsl_interp_accel* left_acc = gsl_interp_accel_alloc();
    gsl_interp_accel* right_acc = gsl_interp_accel_alloc();

    gsl_spline* left_spline = gsl_spline_alloc(gsl_interp_steffen, left_profile.size());
    gsl_spline_init(left_spline, left_profile.data(), left_rv.data(), left_profile.size());

    gsl_spline* right_spline = gsl_spline_alloc(gsl_interp_steffen, right_profile.size());
    gsl_spline_init(right_spline, right_profile.data(), right_rv.data(), right_profile.size());

    for (int i = 0; i < bisector.size(); i++ ) {
        bisector[i] = gsl_spline_eval(left_spline, profile_interp_values[i], left_acc);
    }

    for (int i = 0; i < bisector.size(); i++ ) {
        bisector[i] += gsl_spline_eval(right_spline, profile_interp_values[i], right_acc);
        bisector[i] /= 2.0;
    }

    gsl_spline_free(left_spline);
    gsl_spline_free(right_spline);
    gsl_interp_accel_free(left_acc);
    gsl_interp_accel_free(right_acc);

    return bisector;
}
