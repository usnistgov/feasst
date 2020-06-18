***************************
pyFEASST
***************************

Python support for FEASST includes the following:

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

depend.py
=========

This script is only used for new documentation or plugins.
It generates feasst.i and [plugin]/doc/\*rst files.
