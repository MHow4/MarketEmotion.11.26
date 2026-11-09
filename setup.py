from setuptools import setup, Extension
from Cython.Build import cythonize

ext = Extension(
    "LL_cython",
    ["LL_cython.pyx"],
    extra_compile_args=["-O3", "-fopenmp"],
    extra_link_args=["-fopenmp"],
)
setup(ext_modules=cythonize(ext, annotate=True))