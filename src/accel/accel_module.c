#include "accel_module.h"
#include <math.h>


static PyObject* accel_module_posteriors(PyObject* self, PyObject* args) {
    //ref_base, reads, quals, amp_p_mat, p_ado, f_0
    // OR perhaps have the quals as a 1d np array, the api seems easier
    PyObject *arg0  = NULL,  *arg1  = NULL,  *arg2      = NULL;
    PyObject *reads = NULL,  *quals = NULL,  *amp_p_mat = NULL;
    const char ref_base = 9;
    const double p_ado = 0, f_0 = 0;

    if(!PyArg_ParseTuple(args, "BOOOdd", &ref_base, &arg0, &arg1, &arg2, &p_ado, &f_0)) return NULL;

    reads     = PyArray_FROM_OT(arg0, NPY_BYTE);
    quals     = PyArray_FROM_OT(arg0, NPY_DOUBLE);
    amp_p_mat = PyArray_FROM_OT(arg0, NPY_DOUBLE);

    //When multiple cells held as 2-D array, change to **
    char* read_data = (char*)PyArray_DATA(reads);
    npy_intp n_reads = *PyArray_DIMS(reads);

    //TODO manage threads here

    Py_DECREF(reads);
    Py_DECREF(quals);
    Py_DECREF(amp_p_mat);
    Py_RETURN_NONE;
}


double HWE_prior(int g, double f_0) {
    //Binomial in log space
    double prior = log(f_0) * g;
    prior += (2-g) * log(1-f_0);
    if (g == 1) {
        prior += log(2);
    }
    return prior;
}

static PyMethodDef accel_module_methods[] = { 
    {   
        "posteriors",
        accel_module_posteriors,
        METH_VARARGS,
        "Calculate posterior probabilities in parallel from a cell locus."
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
