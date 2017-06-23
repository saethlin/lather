from setuptools import setup, Extension

install_requires = ['numpy']

lather = Extension('lather',
                   sources=['python_interface.cpp', 'simulation.cpp', 'star.cpp', 'spot.cpp', 'profile.cpp',
                            'inih/ini.c', 'inih/cpp/INIReader.cpp', 'fitrv.cpp', 'boundingshape.cpp',
                            'point.cpp', 'compute_bisector.cpp'],
                   libraries=['gsl', 'gslcblas'],
                   extra_compile_args=['-std=c++14', '-Wno-sign-compare', '-Wno-write-strings'])

setup(name='lather',
      version='0.0.1',
      description='Starspot modeling framework',
      ext_modules=[lather],
      install_requires=install_requires)
