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
#include "fitsim.hpp"


// Fit a simulation to some time and flux data

struct sim_params {
    Simulation sim;
    std::vector<double> time;
    std::vector<double> flux;
    std::vector<bool> spotsArePlage;
};


int sim_f(const gsl_vector* v, void *params, gsl_vector* f) {
    // v is the parameters I'm given
    // params is the struct I'm allowed to pass around
    // f is where I put the new fit values

    // Cast everything out of the struct
    Simulation sim = static_cast<sim_params*>(params)->sim;
    std::vector<double> time = static_cast<sim_params*>(params)->time;
    std::vector<double> flux = static_cast<sim_params*>(params)->flux;
    std::vector<bool> spotsArePlage = static_cast<sim_params*>(params)->spotsArePlage;

    double latitude, longitude, size;

    std::vector<double> output(time.size());
    std::vector<double> dump(time.size());

    // Clear current spots from the simulation
    sim.spots.clear();

    // Add the new spots to the simulation
    for (int i = 0; i < spotsArePlage.size(); i++) {
        latitude = gsl_vector_get(v, 3*i);
        longitude = gsl_vector_get(v, 3*i+1);
        size = gsl_vector_get(v, 3*i+2);
        sim.addSpot(latitude, longitude, size, spotsArePlage[i]);
    }

    sim.observe(time, output, dump, false); // Observe without RV calculation

    // Store residuals in the passed vector
    for (int i = 0; i < time.size(); i++) {
        gsl_vector_set(f, i, output[i]-flux[i]);
    }

    return GSL_SUCCESS;
}


void fit_sim(Simulation* sim, std::vector<double>& time_vector, std::vector<double>& flux_vector) {
    // Set up the solver
    struct sim_params params_struct;
    int iteration = 0, status = 0;
    int npar = sim->spots.size() * 3; // Latitude, longitude, and size for each per spot

    std::vector<bool> spotsArePlage(sim->spots.size());
    for (int i = 0; i < spotsArePlage.size(); i++) {
        spotsArePlage[i] = sim->spots[i].plage;
    }

    gsl_vector *v = gsl_vector_alloc(npar);
    gsl_multifit_fdfsolver *solver;
    gsl_multifit_function_fdf fitFunc;

    // Create the solver
    solver = gsl_multifit_fdfsolver_alloc(gsl_multifit_fdfsolver_lmsder, time_vector.size(), npar);

    // Stick the fit parameters into the struct that GSL needs
    params_struct.time = time_vector;
    params_struct.flux = flux_vector;
    params_struct.spotsArePlage = spotsArePlage; // Save which spots are plages

    for (int i = 0; i < sim->spots.size(); i++) {
        gsl_vector_set(v, 3*i, sim->spots[i].latitude);
        gsl_vector_set(v, 3*i+1, sim->spots[i].longitude);
        gsl_vector_set(v, 3*i+2, sim->spots[i].size);
    }

    // Set up the gsl function
    fitFunc.f = &sim_f;
    fitFunc.df = NULL;
    fitFunc.n = time_vector.size();
    fitFunc.p = npar;
    fitFunc.params = &params_struct;

    // Assign the fitting function to the solver
    gsl_multifit_fdfsolver_set(solver, &fitFunc, v);

    do {
        iteration++;
        status = gsl_multifit_fdfsolver_iterate(solver);
        status = gsl_multifit_test_delta(solver->dx, solver->x, 0.0, 1e-5);
    } while (status == GSL_CONTINUE && iteration < 10000);

    gsl_vector_free(v);
    gsl_multifit_fdfsolver_free(solver);
}
