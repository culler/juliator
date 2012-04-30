from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import sys, os, glob

src = ['iterator.pyx']

setup(
    name = 'iterator',
    version = '1.0',
    description = 'Python extension for iterating rational maps',
    author = 'Marc Culler',
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension('iterator', sources=src)]
    )
