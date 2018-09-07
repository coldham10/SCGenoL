from distutils.core import setup, Extension

accel_module = Extension('accel_module',
                           sources = ['accel_module.c'])

setup(name='genL Accelerator Module',
      version='1.0',
      description='C module for Python to allow multithreading and efficient implementation',
      ext_modules=[accel_module],
      author='Christopher Oldham')
