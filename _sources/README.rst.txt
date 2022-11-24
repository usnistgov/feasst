*************************
README
*************************

The Free Energy and Advanced Sampling Simulation Toolkit (FEASST) is a free,
open-source, public domain software to conduct molecular and particle-based
simulations with Monte Carlo methods.

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

Monte Carlo simulation techniques

* Metropolis
* Wang-Landau
* Transition-matrix
* Mayer-sampling

Thermodynamic ensembles

* Microcanonical ensemble
* Canonical ensemble
* Grand canonical ensemble
* Temperature and growth expanded ensembles

Monte Carlo trials

* Translation, rotation, crankshaft, pivot
* Rigid cluster rotation and translation
* Configurational bias transfers and partial regrowth
* Dual-cut configurational bias
* Aggregation volume bias
* Reptation
* Branches

Interaction potentials

* Hard spheres
* Lennard-Jones with LRC, cut and force shift
* Patchy particles
* Yukawa and charged interactions
* Ewald summation and 2D slab correction
* Bonds, angles and dihedrals
* TraPPE small molecules and n-alkanes
* Slab, cylindrical and spherical confinement
* Cell list and neighbor list

Modern software

* Interface with C++ or as a Python module
* OpenMP parallelization and prefetching
* Checkpoint files to save, restart and analyze simulations

How to install (e.g., compile the executables).
===============================================

.. code-block:: bash

    [apt/yum/brew] install g++ cmake git
    git clone https://github.com/usnistgov/feasst.git
    mkdir feasst/build
    cd feasst/build
    cmake ..
    make install -j4
    # optional python packages for feasst tutorials
    pip install ../pyfeasst jupyter matplotlib pandas scipy

The executables `fst`, which is used to start a simulation, and `rst`, which is used to restart a simulation, should now be located in `/path/to/feasst/build/bin/`.

Troubleshooting install
------------------------

Please :doc:`/CONTACT` us if you run into an issue not listed below.

CentOS 7
~~~~~~~~~

CMake version is usually too old.
Try the command cmake3 instead of cmake.

Windows 10
~~~~~~~~~~~

* Install Windows subsystem for Linux (Ubuntu 16)
* See Ubuntu 16

Ubuntu 16
~~~~~~~~~~

* Update to CMake 3 (https://cmake.org/download/)

macOS Mojave
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* for omp, try brew install libomp

Cray (NERSC CORI)
~~~~~~~~~~~~~~~~~~

* OpenMP functions apparently do not work unless the cray programming environment is disabled.

Ubuntu 18, 20, 22
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* We are not aware of any issues with these OS.

Documentation for a specific version of FEASST
===============================================

You can access the documentation of a specific version of FEASST as follows.

.. code-block:: bash

    git clone https://github.com/usnistgov/feasst.git
    cd feasst
    git checkout nist-pages
    git log
    # find the commit of your version from git log
    # (e.g., 0.19.0 is a50b4fe943832f012373f23658a9497990d70d21)
    git checkout a50b4fe943832f012373f23658a9497990d70d21
    google-chrome index.html

.. include:: CONTACT.rst

.. include:: DISCLAIMER.rst

.. include:: LICENSE.rst
