#ifndef FITRV_HPP
#define FITRV_HPP


#include <vector>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>


struct gauss_params {
    const double* x;
    const double* y;
    const size_t size;
};


int gauss_f(const gsl_vector* v, void* params, gsl_vector* f);


int gauss_df(const gsl_vector* v, void* params, gsl_matrix* J);


std::vector<double> fit_rv(const std::vector<double>& rv_vector,
                           const std::vector<double>& ccf_vector,
                           const std::vector<double>& guess);


#endif