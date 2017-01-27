import os
from setuptools import setup, Extension

os.environ['CC'] = 'clang++'

lather = Extension('lather',
                   sources=['python_interface.cpp', 'simulation.cpp', 'star.cpp', 'spot.cpp', 'profile.cpp', 'inih/ini.c', 'inih/cpp/INIReader.cpp', 'fitrv.cpp', 'fitsim.cpp'],
                   include_dirs=['/usr/local/include', '/usr/local/include/gsl'],
                   library_dirs=['/usr/local/lib'],
                   libraries=['gsl', 'gslcblas'],
                   extra_compile_args=['-std=c++11', '-O3', '-fno-omit-frame-pointer'])

setup(name='lather',
      version='0.0.1',
      description='A starspot modeling program',
      ext_modules=[lather])
