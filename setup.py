from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import sys, os, glob

src = ['iterator.pyx']

setup(
    name = 'juliator',
    version = '1.0',
    description = 'Multiprocessing Julia set browser.',
    author = 'Marc Culler',
    cmdclass = {'build_ext': build_ext},
    packages=['juliator'],
#    package_data={'juliator': doc_files},
    entry_points = {'console_scripts': ['juliator = juliator.app:main']},
    ext_modules = [Extension('iterator', sources=src)]
    )
