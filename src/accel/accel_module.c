#include <Python.h>
#include <numpy/arrayobject.h>

static PyObject* accel_module_posteriors(PyObject *self, PyObject *args) {
    PyObject* arg0 = NULL;
    PyObject* array = NULL;
    if(!PyArg_ParseTuple(args, "O", &arg0)) return NULL;
    array = PyArray_FROM_OT(arg0, NPY_DOUBLE);
    int ndims = PyArray_NDIM(array);
    printf("%d",ndims);
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
