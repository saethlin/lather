#include "fitrv.hpp"


int gauss_f(const gsl_vector* v, void* params, gsl_vector* f) {
    auto height = gsl_vector_get(v, 0);
    auto centroid = gsl_vector_get(v, 1);
    auto width = gsl_vector_get(v, 2);
    auto offset = gsl_vector_get(v, 3);

    const auto x = static_cast<gauss_params*>(params)->x;
    const auto y = static_cast<gauss_params*>(params)->y;
    const auto npts = static_cast<gauss_params*>(params)->size;

    for (auto i = 0; i < npts; i++) {
        auto fit = height * exp(-(x[i]-centroid)*(x[i]-centroid)/(2*width*width)) + offset;
        gsl_vector_set(f, i, fit-y[i]);
    }

    return GSL_SUCCESS;
}


int gauss_df(const gsl_vector* v, void* params, gsl_matrix* J) {
    auto height = gsl_vector_get(v, 0);
    auto centroid = gsl_vector_get(v, 1);
    auto width = gsl_vector_get(v, 2);

    const auto x = static_cast<gauss_params*>(params)->x;
    const auto npts = static_cast<gauss_params*>(params)->size;

    for (auto i = 0; i < npts; i++) {
        auto fit = height * exp(-(x[i]-centroid)*(x[i]-centroid)/(2*width*width));
        gsl_matrix_set (J, i, 0, fit/height);
        gsl_matrix_set (J, i, 1, fit * (x[i]-centroid)/(width*width));
        gsl_matrix_set (J, i, 2, fit * (x[i]-centroid)*(x[i]-centroid)/(width*width*width));
        gsl_matrix_set (J, i, 3, 1);
    }

    return GSL_SUCCESS;
}


std::vector<double> fit_rv(const std::vector<double>& rv_vector, const std::vector<double>& ccf_vector, std::vector<double>& ansatz) {

    std::vector<double> fit_result(ansatz.size());

    // Stick the fit parameters into the struct that GSL needs
    gauss_params input_data {rv_vector.data(), ccf_vector.data(), rv_vector.size()};

    auto param_vector = gsl_vector_alloc(ansatz.size());
    param_vector->data = ansatz.data();

    // Set up the gsl function
    gsl_multifit_function_fdf fitFunc;
    fitFunc.f = &gauss_f;
    fitFunc.df = &gauss_df;
    fitFunc.n = rv_vector.size();
    fitFunc.p = ansatz.size();
    fitFunc.params = &input_data;

    // Create the solver
    auto solver = gsl_multifit_fdfsolver_alloc(gsl_multifit_fdfsolver_lmsder, rv_vector.size(), ansatz.size());

    // Assign the fitting function to the solver
    gsl_multifit_fdfsolver_set(solver, &fitFunc, param_vector);

    int iteration = 0, status = 0;
    do {
        iteration++;
        gsl_multifit_fdfsolver_iterate(solver);
        status = gsl_multifit_test_delta(solver->dx, solver->x, 0.0, 1e-5);
    } while (status == GSL_CONTINUE && iteration < 10000);

    for (auto i = 0; i < ansatz.size(); i++) {
        fit_result[i] = gsl_vector_get(solver->x, i);
    }

    gsl_vector_free(param_vector);
    gsl_multifit_fdfsolver_free(solver);

    return fit_result;
}
