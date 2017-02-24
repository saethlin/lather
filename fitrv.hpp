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
#include <iostream>


struct gauss_params {
    const double* x;
    const double* y;
    const size_t size;
};


int gauss_f(const gsl_vector* v, void* params, gsl_vector* f);


int gauss_df(const gsl_vector* v, void* params, gsl_matrix* J);


std::vector<double> fit_rv(const std::vector<double>& rv_vector, const std::vector<double>& ccf_vector, std::vector<double>& ansatz);


#endif