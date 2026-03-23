from setuptools import setup, Extension
from Cython.Build import cythonize
from numpy import get_include as numpy_get_include
numpy_include_dir = [numpy_get_include()]





ext_modules = [Extension("CRISPResso2Align", ["CRISPResso2Align.pyx"], include_dirs=numpy_include_dir, extra_compile_args=['-w','-Ofast'] ),
                       ]

ext_modules = cythonize(ext_modules, language_level="3")

setup(
    name='CRISPResso2Align',
    ext_modules=cythonize(ext_modules),
    zip_safe=False,
)