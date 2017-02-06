#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "simulation.hpp"
#include <exception>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <Python.h>
#include <numpy/arrayobject.h>


typedef struct {
    PyObject_HEAD
    Simulation CppSimulation;
} PySimulation;


static int PySimulation_init(PySimulation* self, PyObject* args, PyObject* kwargs) {
    int gridSize, spotResolution;
    const char* filename;

    static char* kwdlist1[] = {"grid_size", "spot_resolution", NULL};
    static char* kwdlist2[] = {"filename", NULL};

    try {
        PyArg_ParseTupleAndKeywords(args, kwargs, "s", kwdlist2, &filename);
        self->CppSimulation = Simulation(filename);
        return 0;
    }
    catch (std::exception& e) {
        PyErr_Clear();
        PyArg_ParseTupleAndKeywords(args, kwargs, "ii", kwdlist1, &gridSize, &spotResolution);
        self->CppSimulation = Simulation(gridSize, spotResolution);
        return 0;
    }
}


static PyObject* PySimulation_set_star(PySimulation* self, PyObject* args, PyObject* kwargs) {
    double radius, period, inclination, temperature, spotTempDiff, limbLinear, limbQuadratic;

    static char* kwdlist[] = {"radius", "period", "inclination", "temperature", "spot_temp_diff", "limb_linear", "limb_quadratic", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ddddddd", kwdlist, &radius, &period, &inclination, &temperature, &spotTempDiff, &limbLinear, &limbQuadratic)) {
        return NULL;
    }

    // Set the C++ simulation's star
    self->CppSimulation.setStar(radius, period, inclination, temperature, spotTempDiff, limbLinear, limbQuadratic);

    Py_RETURN_NONE;
}


static PyObject* PySimulation_add_spot(PySimulation* self, PyObject* args, PyObject* kwargs) {
    double latitude, longitude, size;
    bool plage;

    static char* kwdlist[] = {"latitude", "longitude", "size", "plage", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "dddp", kwdlist, &latitude, &longitude, &size, &plage)) {
        return NULL;
    }

    // Add the spot to the C++ simulation
    self->CppSimulation.addSpot(latitude, longitude, size, plage);

    Py_RETURN_NONE;
}


static PyObject* PySimulation_clear_spots(PySimulation* self) {
    self->CppSimulation.clear_spots();
    Py_RETURN_NONE;
}


static PyObject* PySimulation_observe(PySimulation* self, PyObject *args, PyObject *kwargs) {
    PyArrayObject* timeArg = NULL;
    double wavelength;
    bool observeRV;
    static char* kwdlist[] = {"time", "wavelength", "observe_rv", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!d|p", kwdlist,
                                     &PyArray_Type, &timeArg, &wavelength, &observeRV)) {
        return NULL;
    }

    // Create an std::vector that points to the same data as the input array
    npy_intp* dims = PyArray_DIMS(timeArg);
    auto data = (double*)PyArray_DATA(timeArg);
    std::vector<double> time(data, data+dims[0]);

    std::vector<double> flux(time.size());
    std::vector<double> rv(time.size());

    self->CppSimulation.observe(time, flux, rv, wavelength, observeRV);


    /*
    auto dispersion = inplace_phase_fold(data, periods);

    PyObject* output = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    auto output_data = (double*)PyArray_DATA((PyArrayObject*)output);
    for (auto i = 0; i < dispersion.size(); i++) {
        output_data[i] = dispersion[i];
    }
*/

    // Copy std::vector outputs into a dict of numpy arrays
    PyObject* output_flux = PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    PyObject* output_rv = PyArray_SimpleNew(1, dims, NPY_DOUBLE);

    double* output_flux_data = (double*)PyArray_DATA((PyArrayObject*)output_flux);
    double* output_rv_data = (double*)PyArray_DATA((PyArrayObject*)output_rv);
    for (int i = 0; i < dims[0]; i++) {
        output_rv_data[i] = rv[i];
        output_flux_data[i] = flux[i];
    }

    // Build a dictionary to hold the ndarrays
    PyObject* results = PyDict_New();

    PyObject* flux_key = Py_BuildValue("s", "flux");
    PyDict_SetItem(results, flux_key, output_flux);

    PyObject* rv_key = Py_BuildValue("s", "rv");
    PyDict_SetItem(results, rv_key, output_rv);

    return results;
}

/*
static PyObject* PySimulation_fit(PySimulation* self, PyObject *args, PyObject* kwargs) {
    PyObject* timeArg = NULL;
    PyObject* fluxArg = NULL;
    static char* kwdlist[] = {"time", "flux", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO", kwdlist, &timeArg, &fluxArg)) {
        return NULL;
    }

    PyObject* timeArr = PyArray_FROM_OTF(timeArg, NPY_DOUBLE, NPY_IN_ARRAY);
    if (timeArr == NULL) {
        return NULL;
    }

    PyObject* fluxArr = PyArray_FROM_OTF(fluxArg, NPY_DOUBLE, NPY_IN_ARRAY);
    if (fluxArr == NULL) {
        return NULL;
    }

    // Create an std::vector that points to the same data as the input array
    npy_intp* dims = PyArray_DIMS(timeArr);
    double* timeData = (double*)PyArray_DATA(timeArr);
    std::vector<double> time(timeData, timeData+dims[0]);

    double* fluxData = (double*)PyArray_DATA(fluxArr);
    std::vector<double> flux(fluxData, fluxData+dims[0]);

    self->CppSimulation.fit(time, flux);
    Py_RETURN_NONE;
}
*/

/*
static PyObject* PySimulation_str(PySimulation* self) {
    std::ostringstream reprStream;
    reprStream << std::boolalpha;

    reprStream << "Simulation:\n";
    reprStream << "    Grid size: " << self->CppSimulation.gridSize << "\n";
    reprStream << "    Spot resolution: " << self->CppSimulation.spotResolution << "\n";
    reprStream << "Star:\n";
    reprStream << "    Radius: " << self->CppSimulation.star.radius << "\n";
    reprStream << "    Period: " << self->CppSimulation.star.period << "\n";
    reprStream << "    Inclination: " << self->CppSimulation.star.inclination << "\n";
    reprStream << "    Temperature: " << self->CppSimulation.star.temperature << "\n";
    reprStream << "    Spot temperature difference: " << self->CppSimulation.star.spotTempDiff << "\n";
    reprStream << "    Limb darkening coefficients: " << self->CppSimulation.star.limbLinear << " " << self->CppSimulation.star.limbQuadratic << "\n";

    std::vector<Spot>::iterator spot;
    for (spot = self->CppSimulation.spots.begin(); spot != self->CppSimulation.spots.end(); ++spot) {
        reprStream << "Spot:\n";
        reprStream << "    Latitude: " << spot->latitude << "\n";
        reprStream << "    Longitude: " << spot->longitude << "\n";
        reprStream << "    Size: " << spot->size << "\n";
        reprStream << "    Plage: " << spot->plage << "\n";
    }

    return PyUnicode_FromString(reprStream.str().c_str());
}
*/


static PyMethodDef PySimulation_methods[] = {
    {"set_star", (PyCFunction)PySimulation_set_star, METH_KEYWORDS|METH_VARARGS,
     "Set the star parameters"},
    {"add_spot", (PyCFunction)PySimulation_add_spot, METH_KEYWORDS|METH_VARARGS,
     "Add a spot to the simulation"},
    {"clear_spots", (PyCFunction)PySimulation_clear_spots, METH_NOARGS,
     "Remove all spots on the simulation"},
    {"observe", (PyCFunction)PySimulation_observe, METH_KEYWORDS|METH_VARARGS,
     "Compute simulated observations at the given times"},
    //{"fit", (PyCFunction)PySimulation_fit, METH_KEYWORDS|METH_VARARGS,
    // "Attempt to fit data with the current simulation"},
    //{"toString", (PyCFunction)PySimulation_str, METH_NOARGS,
    // "Produce a string representation"},
    {NULL}  /* Sentinel */
};


static PyTypeObject PySimulationType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "lather.Simulation",       /* tp_name */
    sizeof(PySimulation),      /* tp_basicsize */
    0,                         /* tp_itemsize */
    0,                         /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_reserved */
    0, //(reprfunc)PySimulation_str,          /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash  */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT |
        Py_TPFLAGS_BASETYPE,   /* tp_flags */
    "lather Simulation objects",/* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    PySimulation_methods,      /* tp_methods */
    0, //Simulation_members,   /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)PySimulation_init,/* tp_init */
    0,                         /* tp_alloc */
    PyType_GenericNew,         /* tp_new */
};


static PyMethodDef latherMethods[] = {
    {NULL, NULL, 0, NULL}        // Sentinel
};


static struct PyModuleDef latherModule = {
   PyModuleDef_HEAD_INIT,
   "lather", // name of module
   NULL, // module documentation, may be NULL
   -1,  //size of per-interpreter state of the module, or -1 if the module keeps state in global variables
   latherMethods
};


PyMODINIT_FUNC PyInit_lather(void) {
    import_array();

    PyObject* module;

    PySimulationType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&PySimulationType) < 0) {
        return NULL;
    }

    module = PyModule_Create(&latherModule);
    if (module == NULL) {
        return NULL;
    }

    Py_INCREF(&PySimulationType);
    PyModule_AddObject(module, "Simulation", (PyObject*)&PySimulationType);

    return module;
}