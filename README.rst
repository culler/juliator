Juliator
========

Description
-----------

Juliator is a graphical application for viewing Mandelbrot and
Julia sets for the usual quadratic dynamical system z -> z^2 + a.
The program is written in Python using tkinter and PIL for the
graphics.  The computations are done in a C extension module
wrapped using Cython.  An unusual feature is that the computation
is done in parallel with a separate process for each CPU core
and, more interesingly, these processes use shared memory rather
than communicating via an IPC scheme.  Each process is assigned
to draw a section of a common bitmap image.

Installation
------------

For macOS there is a disk image available in the Downloads
section which contains a standalone app suitable for copying
to the applications folder.

Other platforms should install with

::
  python3 setup.py install

or

::
  sudo python3 setup.py install

and then run the app with

::
  python3 -m juliator.app

  