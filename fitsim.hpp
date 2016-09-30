#ifndef FITSIM_HPP
#define FITSIM_HPP
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_vector.h>
#include "simulation.hpp"

void fit_sim(Simulation* sim, std::vector<double>& time, std::vector<double>& flux);

int sim_f(const gsl_vector *v, void *params, gsl_vector *f);

#endif