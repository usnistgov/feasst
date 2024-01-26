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

The features available in the :doc:`/plugin/text_interface` are summarized as follows:

Monte Carlo simulation techniques

* Metropolis
* Wang-Landau
* Transition-matrix
* Mayer-sampling

Thermodynamic ensembles

* Microcanonical ensemble
* Canonical ensemble
* Grand-canonical ensemble
* Temperature and growth expanded ensembles
* Gibbs ensemble

Interaction potentials

* Lennard-Jones and Mie with cut/force shift and corrections
* Hard spheres and square wells
* Patchy and anisotropic particles
* Yukawa and charged interactions
* Ewald summation and 2D slab correction
* Bonds, angles and dihedrals
* TraPPE small molecules and n-alkanes
* Slab, cylindrical and spherical confinement
* Cell list and neighbor list

Monte Carlo trials

* Translation, rotation, crankshaft and pivot
* Rigid cluster rotation and translation
* Configurational bias transfers and partial regrowth
* Dual-cut configurational bias
* Aggregation volume bias
* Reptation
* Branches

Modern software

* Interface as text input, C++ or Python module
* Server interface to C++ or Python clients
* OpenMP parallelization and prefetching
* Checkpoint files to save, restart and analyze simulations

How to get started
===============================================

First, compile the executables as described below.
Second, find :doc:`tutorial/README` that are closest to what you would like to accomplish.
Third, reproduce the expected result of those tutorials.
Fourth, use the :doc:`../plugin/text_interface` documentation to better understand and modify the tutorial to accomplish your goals.
When you :doc:`../CONTACT` us with issues, include the text input file instead of the python script that generates that text file.

How to install (i.e., compile the executables)
===============================================

.. code-block:: bash

    [apt/yum/brew] install g++ cmake git python3
    git clone https://github.com/usnistgov/feasst.git
    mkdir feasst/build
    cd feasst/build
    cmake ..
    make install -j4
    # optional python packages for feasst tutorials
    pip install ../pyfeasst jupyter matplotlib pandas scipy

The executables `fst` and `rst` should appear in `/path/to/feasst/build/bin/`.
Text input files are run using `fst < input.txt` while simulations are restarted using `rst checkpoint.txt`.
It is important to provide pip a path to the specific pyfeasst directory in feasst to ensure the versions match (e.g., do not leave out the "../" above).

Troubleshooting install
------------------------

Please :doc:`/CONTACT` us if you run into an issue not listed below.

CentOS 7
~~~~~~~~~

CMake version is usually too old.
Try the command cmake3 instead of cmake.

Rocky 8
~~~~~~~~

* yum install gcc-c++

Ubuntu 16
~~~~~~~~~~

* Update to CMake 3 (https://cmake.org/download/)

Ubuntu 18, 20, 22
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* We are not aware of any install issues with these OS.

Cray (NERSC CORI)
~~~~~~~~~~~~~~~~~~

* OpenMP functions may not work unless the cray programming environment is disabled.

macOS Mojave
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* for omp, try brew install libomp

Windows 10
~~~~~~~~~~~

* Install Windows subsystem for Linux (Ubuntu 16)
* See Ubuntu 16

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
