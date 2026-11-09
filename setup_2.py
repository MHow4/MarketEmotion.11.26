import numpy as np
from setuptools import setup, Extension
from Cython.Build import cythonize

ext = Extension(
    "LL_cython",
    ["LL_cython.pyx"],
    include_dirs=[np.get_include()],
    extra_compile_args=["/O2", "/openmp"],
    extra_link_args=["/openmp"],
)
setup(ext_modules=cythonize(ext, annotate=True))