#include <Python.h>
#include <numpy/arrayobject.h>

static PyObject* accel_module_posteriors(PyObject *self, PyObject *args)
{
    printf("Hello Worldd\n");
    Py_RETURN_NONE;
}

static PyMethodDef accel_module_methods[] = { 
    {   
        "posteriors",
        accel_module_posteriors,
        METH_VARARGS,
        "Print 'hello world' from a method defined in a C extension."
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

PyMODINIT_FUNC PyInit_accel_module(void)
{
    Py_Initialize();
    import_array();

    return PyModule_Create(&accel_module_definition);
}
