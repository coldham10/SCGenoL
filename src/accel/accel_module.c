#include "accel_module.h"
#include <math.h>

static double HWE_prior(int g, double f_0) {
    //Binomial in log space
    double prior = log(f_0) * g;
    prior += (2-g) * log(1-f_0);
    if (g == 1) {
        prior += log(2);
    }
    return prior;
}

static PyObject* accel_module_posteriors(PyObject* self, PyObject* args) {
    PyObject* arg0 = NULL;
    PyObject* amp_p_mat = NULL;
    if(!PyArg_ParseTuple(args, "O", &arg0)) return NULL;
    amp_p_mat = PyArray_FROM_OT(arg0, NPY_DOUBLE);
    int ndims = PyArray_NDIM(amp_p_mat);
    printf("%d\n", ndims);
    double test_prior = HWE_prior(1, 0.5);
    printf("%f\n", test_prior);
    Py_DECREF(amp_p_mat);
    Py_RETURN_NONE;
}


static PyMethodDef accel_module_methods[] = { 
    {   
        "posteriors",
        accel_module_posteriors,
        METH_VARARGS,
        "Calculate posterior probabilities in parallel from a group of cell loci"
    },  
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef accel_module_definition = { 
    PyModuleDef_HEAD_INIT,
    "accel_module",
    "A module to allow eficient multiprocessing for computationaly intensive tasks in genL",
    -1, 
    accel_module_methods
};

PyMODINIT_FUNC PyInit_accel_module(void) {
    Py_Initialize();
    import_array();

    return PyModule_Create(&accel_module_definition);
}
