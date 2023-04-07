from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

ext_modules=[
    Extension("supercoiling_cython",
              sources=["supercoiling_cython.pyx"],
              libraries=["m"] # Unix-like specific
    )
]

setup(
    name = "supercoiling_cython",
    ext_modules = cythonize(ext_modules),
    include_dirs=[numpy.get_include()]
)
