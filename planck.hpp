#ifndef LATHER_PLANCK_H
#define LATHER_PLANCK_H


#include <math.h>
#include <gsl/gsl_integration.h>


const double c = 299792458;
const double h = 6.62606896e-34;
const double k_b = 1.380e-23;


inline double planck(double wavelength, void* params) {
    double temperature = *(double*)params;
    return 2*h*c*c*1./(wavelength*wavelength*wavelength*wavelength*wavelength*(expm1((h*c)/(wavelength*k_b*temperature))));
}


inline double planck_integral(double temperature, double wave_min, double wave_max) {
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

    double result, error;

    gsl_function F;
    F.function = &planck;
    F.params = &temperature;

    gsl_integration_qags(&F, wave_min, wave_max, 0, 1e-7, 1000,
                         w, &result, &error);

    gsl_integration_workspace_free(w);

    return result;
}

#endif