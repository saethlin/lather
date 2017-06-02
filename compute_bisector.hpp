#ifndef LATHER_COMPUTE_BISECTOR_HPP
#define LATHER_COMPUTE_BISECTOR_HPP

#include <vector>
#include <algorithm>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <iostream>

std::vector<double> compute_bisector(const std::vector<double>& rv, const std::vector<double>& profile);

#endif //LATHER_COMPUTE_BISECTOR_HPP
