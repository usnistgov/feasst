*************************
README
*************************

.. image:: https://travis-ci.com/hhatch/feasst.svg?branch=master
    :target: https://travis-ci.com/hhatch/feasst

.. image:: https://codecov.io/gh/hhatch/feasst/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/hhatch/feasst

.. image:: https://www.codefactor.io/repository/github/hhatch/feasst/badge/master
    :target: https://www.codefactor.io/repository/github/hhatch/feasst/overview/master
    :alt: CodeFactor

The Free Energy and Advanced Sampling Simulation Toolkit (FEASST) is a free,
open-source, modular program to conduct molecular and particle-based
simulations with flat-histogram Monte Carlo and molecular dynamics methods.

.. note::

   Manuscript: https://doi.org/10.6028/jres.123.004

   Website: https://pages.nist.gov/feasst/

   Website DOI: https://doi.org/10.18434/M3S095

   Code repository: https://github.com/usnistgov/feasst

   Discussion list: https://groups.google.com/a/list.nist.gov/d/forum/feasst

Features
================================================================================

.. image:: dev/sphinx/feasst.png
   :target: https://pages.nist.gov/feasst
   :align: right

:doc:`/plugin/README` contains the list of features available.
Some features include but are not limited to the following:

Simulation techniques

* Wang-Landau Monte Carlo
* Transition-matrix Monte Carlo
* Metropolis Monte Carlo
* Mayer-sampling Monte Carlo
* Configurational bias

Thermodynamic ensembles

* Microcanonical ensemble
* Canonical ensemble
* Grad canonical ensemble

Intermolecular interactions

* Hard spheres
* Lennard-Jones with Yukawa, LRC, force shift
* Patchy particles
* Charged interactions with the Ewald summation
* Cylindrical and slit pore confinement

Modern software

* Interface with C++ or as a Python module
* OpenMP parallelization
* Checkpointing to save and restart simulations

Usage
================================================================================

The following example Lennard-Jones Monte Carlo simulation may be found in :doc:`/tutorial/README`.

.. literalinclude:: tutorial/lj_brief.py
   :language: py

Python usage
================================================================================

.. code-block:: bash

    mkdir /path/to/feasst/build
    cd /path/to/feasst/build
    cmake -DUSE_SWIG=ON ..
    make _feasst -j12
    ../py/run.sh ../py/test.py  # optional test

CMake attempts to find the python libraries.
But you may want to specify them manually as follows:

.. code-block:: bash

    cmake -DUSE_SWIG=ON -DSET_PYTHON_PATH=ON -DPYTHON_INCLUDE_PATH=/path/to/anaconda/include/python3.6m -DPYTHON_LIBRARIES=/path/to/anaconda/lib/libpython3.6m.so ..

It is recommended to use SWIG version 3 or 4 and not 2.

Python scrips are then run using /path/to/feasst/py/run.sh, which sets PYTHONPATH to /path/to/feasst/build.
Alternatively, you could set the python path manually.

C++ usage
================================================================================

First, install the C++ library.

.. code-block:: bash

    mkdir /path/to/feasst/build
    cd /path/to/feasst/build
    cmake ..
    make install -j12
    make test         # optional test

Then, compile the specific simulation you wish to run (e.g., tutorial).

.. code-block:: bash

    cd /path/to/feasst/tutorial/
    mkdir build; cd build
    cmake ..
    make
    ./lj_brief

CMake defaults to install in the build directory.
But you can also specify the path as follows.

.. code-block:: bash

    cmake -DCMAKE_INSTALL_PREFIX=/path/to/install/dir ..

Later when you build your executable you need to specify this path as follows:

.. code-block:: bash

    cmake -DCMAKE_PREFIX_PATH=/path/to/install/dir ..

.. include:: CONTACT.rst

.. include:: DISCLAIMER.rst

.. include:: LICENSE.rst
