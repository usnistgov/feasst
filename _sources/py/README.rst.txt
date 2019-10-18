***************************
pyFEASST
***************************

Python support for FEASST includes the following:

run.sh
=========

Shell script which assumes that feasst was compiled in the directory /path/to/feasst/build/ and that run.sh itself is in /path/to/feasst/py.
The script sets PYTHONPATH.
This script assumes python3 but you could easily modify it for your use.
CMakeLists.txt also assumes python3 for finding LIBS, so you'd have to manually specify the python paths as described in the installation instructions.

pyfeasst.py
============

Collection of utility python functions and classes for use with FEASST.

test.py
========

Unittests for python interface.

feasst.i
=========

This is the SWIG interface file for wrapping C++ to Python.
Typically only developers would be interested in this file.
It is generated using a script because the order of including headers is important.

