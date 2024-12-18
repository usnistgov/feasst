*************************
Overview
*************************

.. include:: README.rst

.. include:: build/build.rst

How to learn more
===============================================

* Review the :doc:`first <tutorial/tutorial>` and :doc:`second tutorial<tutorial/launch>`.

  * Copy/paste or use the URL to find the code (e.g., https://pages.nist.gov/feasst/tutorial/launch.html is ``$HOME/feasst/tutorial/launch.py``).
  * See ``python launch.py --help`` (e.g., adjust ``--feasst_install`` or ``--hours_terminate``).
* Find :doc:`tutorial/README` that are closest to what you would like to accomplish.
* Reproduce the expected result of those :doc:`tutorial/README`.
* To modify :doc:`tutorial/README` to accomplish your goals, refer to the :doc:`../plugin/text_interface` documentation.
* Compare the energy of a :doc:`reference configuration<plugin/monte_carlo/tutorial/tutorial_0_ref_configs>` with a trusted source to ensure the model, :doc:`/particle/README` return an expected result.
* :doc:`/CONTACT` us.

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

Troubleshooting install
------------------------

Please :doc:`/CONTACT` us if you run into an issue not listed below.

Ubuntu 18, 20, 22, 24 and Rocky 8 and 9
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* We are not aware of any install issues with these OS.

Ubuntu 16
~~~~~~~~~~

* Update to CMake 3 (https://cmake.org/download/)

CentOS 7
~~~~~~~~~

CMake version is usually too old.
Try the command cmake3 instead of cmake.

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

.. include:: CONTACT.rst

.. include:: DISCLAIMER.rst

.. include:: LICENSE.rst
