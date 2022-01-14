"""
from setuptools import setup
from Cython.Build import cythonize

python setup_infection_model_URT.py build_ext --inplace
"""
from setuptools import setup
from Cython.Build import cythonize
import numpy

setup(
   ext_modules = cythonize("infection_model_cython_URT.pyx"),
   include_dirs=[numpy.get_include()]
)
