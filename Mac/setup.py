import setuptools

# Hack to patch setuptools so that it treats Cython
# as a replacement for pyrex.

from distutils.core import Extension as _Extension
from setuptools.dist import _get_unpatched
_Extension = _get_unpatched(_Extension)
try:
    from Cython.Distutils import build_ext
except ImportError:
    have_cython = False
else:
    have_cython = True

class Extension(_Extension):
    """
    This modified version of setuptools Extension allows us
    to use Cython instead of pyrex.  If Cython is not installed
    on this system, it will assume that a Cython-generated .c
    file is present in the distribution.
    """
    if not have_cython:
        # convert .pyx extensions to .c
        def __init__(self,*args,**kw):
            _Extension.__init__(self,*args,**kw)
            sources = []
            for s in self.sources:
                if s.endswith('.pyx'):
                    sources.append(s[:-3]+'c')
                else:
                    sources.append(s)
            self.sources = sources

class Library(Extension):
    """Just like a regular Extension, but built as a library instead"""

import sys, distutils.core, distutils.extension
distutils.core.Extension = Extension
distutils.extension.Extension = Extension
if 'distutils.command.build_ext' in sys.modules:
    sys.modules['distutils.command.build_ext'].Extension = Extension
# End of hack

from setuptools import setup
from pkg_resources import load_entry_point

SRC = ['iterator.pyx']
APP = ['Juliator.py']
OPTIONS = {'argv_emulation': False,
           'iconfile': 'Juliator.icns',
           'includes': 'multiprocessing',
}

setup(
    name = 'juliator',
    version = '1.0',
    description = 'Multiprocessing Julia set browser.',
    author = 'Marc Culler',
    cmdclass = {'build_ext': build_ext},
    app=APP,
    packages=['juliator'],
#    package_data={'juliator': doc_files},
    options={'py2app': OPTIONS},
    entry_points = {'console_scripts': ['juliator = juliator.app:main']},
    ext_modules = [Extension('iterator', sources=SRC)],
    setup_requires=['py2app']
    )
