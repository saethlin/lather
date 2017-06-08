from setuptools import setup, Extension
import os

install_requires = ['numpy']

os.environ['CC'] = 'g++'

lather = Extension('lather',
                   sources=['python_interface.cpp', 'simulation.cpp', 'star.cpp', 'spot.cpp', 'profile.cpp',
                            'inih/ini.c', 'inih/cpp/INIReader.cpp', 'fitrv.cpp', 'boundingshape.cpp',
                            'point.cpp', 'compute_bisector.cpp'],
                   include_dirs=['/usr/local/include', '/usr/local/include/gsl'],
                   library_dirs=['/usr/local/lib'],
                   libraries=['gsl', 'gslcblas'],
                   extra_compile_args=['-std=c++14', '-O3'])

setup(name='lather',
      version='0.0.1',
      description='Starspot modeling framework',
      ext_modules=[lather],
      install_requires=install_requires)
