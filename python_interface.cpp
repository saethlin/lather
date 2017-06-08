#include <exception>
#include <stdio.h>
#include <sstream>
#include <vector>
#include <stdbool.h>
#include <Python.h>
#include <numpy/arrayobject.h>
#include <memory>
#include "simulation.hpp"


typedef struct {
    PyObject_HEAD
    std::unique_ptr<Simulation> CppSimulation;
} PySimulation;


static int PySimulation_init(PySimulation* self, PyObject* args, PyObject* kwargs) {
    const char* filename;
    static char* kwdlist2[] = {"filename", NULL};

    try {
        PyArg_ParseTupleAndKeywords(args, kwargs, "s", kwdlist2, &filename);
        self->CppSimulation = std::make_unique<Simulation>(filename);
        return 0;
    }
    catch (std::exception& e) {
        PyErr_Clear();
        self->CppSimulation = std::make_unique<Simulation>();
        return 0;
    }
}


static PyObject* PySimulation_set_star(PySimulation* self, PyObject* args, PyObject* kwargs) {
    int grid_size;
    double radius, period, inclination, temperature, spotTempDiff, limbLinear, limbQuadratic;

    static char* kwdlist[] = {"grid_size", "radius", "period", "inclination", "temperature", "spot_temp_diff", "limb_linear", "limb_quadratic", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "iddddddd", kwdlist, &grid_size, &radius, &period, &inclination, &temperature, &spotTempDiff, &limbLinear, &limbQuadratic)) {
        return NULL;
    }

    // Set the C++ simulation's star
    self->CppSimulation->set_star(grid_size, radius, period, inclination, temperature, spotTempDiff, limbLinear, limbQuadratic);

    Py_RETURN_NONE;
}


static PyObject* PySimulation_add_spot(PySimulation* self, PyObject* args, PyObject* kwargs) {
    double latitude, longitude, size;
    bool plage = false;
    bool mortal = false;

    static char* kwdlist[] = {"latitude", "longitude", "size", "plage", "mortal", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "dddpp", kwdlist, &latitude, &longitude, &size, &plage, &mortal)) {
        return NULL;
    }

    self->CppSimulation->add_spot(latitude, longitude, size, plage, mortal);

    Py_RETURN_NONE;
}


static PyObject* PySimulation_clear_spots(PySimulation* self) {
    self->CppSimulation->clear_spots();
    Py_RETURN_NONE;
}


static PyObject* PySimulation_observe_rv(PySimulation* self, PyObject *args, PyObject *kwargs) {
    PyObject* timeArg;
    double wavelength_min, wavelength_max;
    static char* kwdlist[] = {"time", "wavelength_min", "wavelength_max", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!dd", kwdlist,
                                     &PyArray_Type, &timeArg,
                                     &wavelength_min, &wavelength_max)) {
        return NULL;
    }

    // Create an std::vector that points to the same data as the input array
    npy_intp* dims = PyArray_DIMS(timeArg);
    double* data = (double*)PyArray_DATA(timeArg);
    std::vector<double> time(data, data+dims[0]);

    auto rv_observations = self->CppSimulation->observe_rv(time, wavelength_min, wavelength_max);

    // Construct the output array for the rv values
    PyObject* output_rv = PyArray_SimpleNew(1, dims, NPY_DOUBLE);

    double* output_rv_data = (double*)PyArray_DATA(output_rv);
    for (int i = 0; i < dims[0]; i++) {
        output_rv_data[i] = rv_observations[i].rv;
    }

    // Construct the output array for the bisectors
    npy_intp bis_dims[2] = {time.size(), rv_observations[0].bisector.size()};
    PyObject* output_bisectors = PyArray_SimpleNew(2, bis_dims, NPY_DOUBLE);

    double* bisector_data = (double*)PyArray_DATA(output_bisectors);
    int i = 0;
    for (const auto& observation : rv_observations) {
        for (const auto& val : observation.bisector) {
            bisector_data[i] = val;
            i++;
        }
    }

    PyObject* output_tuple = Py_BuildValue("OO", output_rv, output_bisectors);

    return output_rv;
}


static PyObject* PySimulation_observe_flux(PySimulation* self, PyObject *args, PyObject *kwargs) {
    PyObject* timeArg = NULL;
    double wavelength_min, wavelength_max;
    static char* kwdlist[] = {"time", "wavelength_min", "wavelength_max", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!dd", kwdlist,
                                     &PyArray_Type, &timeArg,
                                     &wavelength_min, &wavelength_max)) {
        return NULL;
    }

    // Create an std::vector that points to the same data as the input array
    npy_intp* dims = PyArray_DIMS(timeArg);
    double* data = (double*)PyArray_DATA(timeArg);
    std::vector<double> time(data, data+dims[0]);

    auto flux = self->CppSimulation->observe_flux(time, wavelength_min, wavelength_max);

    PyObject* output_flux = PyArray_SimpleNew(1, dims, NPY_DOUBLE);

    double* output_flux_data = (double*)PyArray_DATA(output_flux);
    for (int i = 0; i < dims[0]; i++) {
        output_flux_data[i] = flux[i];
    }
    return output_flux;
}


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
    {"observe_rv", (PyCFunction)PySimulation_observe_rv, METH_KEYWORDS|METH_VARARGS,
     "Compute simulated observations at the given times"},
    {"observe_flux", (PyCFunction)PySimulation_observe_flux, METH_KEYWORDS|METH_VARARGS,
     "Compute simulated observations at the given times"},
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