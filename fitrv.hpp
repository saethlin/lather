#ifndef FITRV_HPP
#define FITRV_HPP
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_vector.h>


void fit_rv(std::vector<double>& rv, std::vector<double>& ccf, std::vector<double>& fit_output);

int gauss_df(const gsl_vector* v, void* p, gsl_matrix* J);

int gauss_f(const gsl_vector *v, void *params, gsl_vector *f);

#endif