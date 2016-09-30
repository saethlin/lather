#include <iostream>
#include <math.h>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_vector.h>
#include "fitrv.hpp"


// Fit a profile with a gaussian to determine the radial velocity

struct gauss_params {
    double* x;
    double* y;
    int npts;
    //std::vector<double> x;
    //std::vector<double> y;
};


int gauss_f(const gsl_vector *v, void *params, gsl_vector *f) {
    double fit;
    double height = gsl_vector_get(v, 0);
    double centroid = gsl_vector_get(v, 1);
    double width = gsl_vector_get(v, 2);
    double offset = gsl_vector_get(v, 3);

    double* x = static_cast<gauss_params*>(params)->x;
    double* y = static_cast<gauss_params*>(params)->y;
    int& npts = static_cast<gauss_params*>(params)->npts;
    //std::vector<double>& x = static_cast<gauss_params*>(params)->x;
    //std::vector<double>& y = static_cast<gauss_params*>(params)->y;

    //for (int i = 0; i < x.size(); i++) {
    for (int i = 0; i < npts; i++) {
        fit = height * exp(-(x[i]-centroid)*(x[i]-centroid)/(2*width*width)) + offset;
        gsl_vector_set(f, i, fit-y[i]);
    }

    return GSL_SUCCESS;
}


int gauss_df(const gsl_vector* v, void* params, gsl_matrix* J) {
    double height = gsl_vector_get(v, 0);
    double centroid = gsl_vector_get(v, 1);
    double width = gsl_vector_get(v, 2);
    double fit;

    double* x = static_cast<gauss_params*>(params)->x;
    int& npts = static_cast<gauss_params*>(params)->npts;
    //std::vector<double>& x = static_cast<gauss_params*>(params)->x;

    //for (int i = 0; i < x.size(); i++) {
    for (int i = 0; i < npts; i++) {
        fit = height * exp(-(x[i]-centroid)*(x[i]-centroid)/(2*width*width));
        gsl_matrix_set (J, i, 0, fit/height);
        gsl_matrix_set (J, i, 1, fit * (x[i]-centroid)/(width*width));
        gsl_matrix_set (J, i, 2, fit * (x[i]-centroid)*(x[i]-centroid)/(width*width*width));
        gsl_matrix_set (J, i, 3, 1);
    }

    return GSL_SUCCESS;
}


void fit_rv(std::vector<double>& rv_vector, std::vector<double>& ccf_vector, std::vector<double>& fit_result) {
    // Set up the solver
    struct gauss_params params_struct;
    int iteration = 0, status = 0, npar = 4;

    // Separate data from the vectors for compatabilty with older versions of gsl
    double* rv = &rv_vector[0];
    double* ccf = &ccf_vector[0];
    unsigned int npts = rv_vector.size();

    gsl_vector *v = gsl_vector_alloc(npar);
    gsl_multifit_fdfsolver *solver;
    gsl_multifit_function_fdf fitFunc;

    // Create the solver
    solver = gsl_multifit_fdfsolver_alloc(gsl_multifit_fdfsolver_lmsder, npts, npar);

    // Stick the fit parameters into the struct that GSL needs
    params_struct.x = rv;
    params_struct.y = ccf;
    params_struct.npts = npts;

    for (int i = 0; i < 4; i++) {
        gsl_vector_set(v, i, fit_result[i]);
    }

    // Set up the gsl function
    fitFunc.f = &gauss_f;
    fitFunc.df = &gauss_df;
    fitFunc.n = npts;
    fitFunc.p = npar;
    fitFunc.params = &params_struct;

    // Assign the fitting function to the solver
    gsl_multifit_fdfsolver_set(solver, &fitFunc, v);

    do {
        iteration++;
        status = gsl_multifit_fdfsolver_iterate(solver);
        status = gsl_multifit_test_delta(solver->dx, solver->x, 0.0, 1e-5);
    } while (status == GSL_CONTINUE && iteration < 10000);

    for (int i = 0; i < 4; i++) {
        fit_result[i] = gsl_vector_get(solver->x, i);
    }

    gsl_vector_free(v);
    gsl_multifit_fdfsolver_free(solver);
}
