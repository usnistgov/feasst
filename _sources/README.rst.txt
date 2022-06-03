*************************
README
*************************

The Free Energy and Advanced Sampling Simulation Toolkit (FEASST) is a free,
open-source, modular program to conduct molecular and particle-based
simulations with flat-histogram Monte Carlo methods.

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

Tutorial
================================================================================

`Tutorial <tutorial/README.html>`_ describes an example Lennard-Jones Monte Carlo simulation.
Also check out tutorials for :doc:`/plugin/README`.
In particular, `MonteCarlo <plugin/monte_carlo/README.html>`_ and `FlatHistogram <plugin/flat_histogram/README.html>`_.

The search box for the html documentation was disabled for security reasons.
But the html is generated entirely from the downloaded code.
Thus, the bash command "grep" is a great option to search for more information on classes and their arguments.
For example, if you would like more information on `RandomMT19937 <plugin/math/doc/RandomMT19937.html>`_ but are not sure where to find it, you could search headers files

.. code-block:: bash

   grep -r --include=*.h RandomMT19937

And find that the class is part of the `Math <plugin/math/README.html>`_ plugin.
Searching the GitHub repository is another option.

How to install (e.g., compile the executables).
===============================================

.. code-block:: bash

    [apt/yum/brew] install g++ cmake git
    git clone https://github.com/usnistgov/feasst.git
    mkdir feasst/build
    cd feasst/build
    cmake ..
    make install -j4
    pip install jupyter matplotlib pandas scipy # optional for tutorials

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

Ubuntu 18, 20 and macOS Mojave
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* We are not aware of any issues with these OS.

Cray (NERSC CORI)
~~~~~~~~~~~~~~~~~~

* OpenMP functions apparently do not work unless the cray programming environment is disabled.

.. include:: CONTACT.rst

.. include:: DISCLAIMER.rst

.. include:: LICENSE.rst
