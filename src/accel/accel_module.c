#include "accel_module.h"
#include <math.h>


static PyObject* accel_module_posteriors(PyObject* self, PyObject* args) {
    //ref_base, read_string, qual_list, amp_p_mat, p_ado, f_0
    //TODO its probably more efficient to read quals as string and convert in C to double
    // OR perhaps have the quals as a 1d np array, the api seems easier
    PyObject *arg0=NULL, *arg1=NULL;
    const char ref_base = '\0';
    const char* read_string;
    PyObject* amp_p_mat = NULL;
    if(!PyArg_ParseTuple(args, "BsOO", &ref_base, &read_string, &arg0, &arg1)) return NULL;
    amp_p_mat = PyArray_FROM_OT(arg1, NPY_DOUBLE);
    int ndims = PyArray_NDIM(amp_p_mat);
    printf("%d\n", ndims);
    printf("%d\n", (int)ref_base);
    printf("%s\n", read_string);
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
