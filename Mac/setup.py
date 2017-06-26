import setuptools

from setuptools import setup
from pkg_resources import load_entry_point

SRC = ['iterator.pyx']
APP = ['Juliator.py']
OPTIONS = {'argv_emulation': False,
           'iconfile': 'Juliator.icns',
           'includes': 'multiprocessing',
}

setup(
    name = 'Juliator',
    version = '1.0',
    description = 'Multiprocessing Julia set browser.',
    author = 'Marc Culler',
    app=APP,
    options={'py2app': OPTIONS},
    entry_points = {'console_scripts': ['juliator = juliator.app:main']},
    setup_requires=['py2app']
    )

# Py2app doesn't put the Tcl scripts in a place where tkinter
# can find them.  So we have to fix up the bundle.
