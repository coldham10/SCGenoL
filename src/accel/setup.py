from distutils.core import setup, Extension
import numpy as np

accel_module = Extension('accel_module',
                           sources = ['accel_module.c'],
                           include_dirs=[np.get_include()])

setup(name='genL Accelerator Module',
      version='1.0',
      description='C module for Python to allow multithreading and efficient implementation',
      ext_modules=[accel_module],
      author='Christopher Oldham')
