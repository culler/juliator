# -*- coding: utf-8 -*-

# Copyright© 2012-2017 by Marc Culler and others.
# This file is part of Juliator.
#
# Juliator is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# Juliator is distributed in the hope that it will be useful
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Juliator.  If not, see <http://www.gnu.org/licenses/>.

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
