*************************
Overview
*************************

The Free Energy and Advanced Sampling Simulation Toolkit (FEASST) is a free and open-source software to conduct molecular- and particle-based simulations with Monte Carlo methods.
New users can start with the `website <https://pages.nist.gov/feasst/>`_ (`DOI <https://doi.org/10.18434/M3S095>`_), `manuscript <https://doi.org/10.1063/5.0224283>`_, `GitHub discussion <https://github.com/usnistgov/feasst/discussions>`_ and a `five minute video <https://www.nist.gov/video/how-use-feasst-0255-monte-carlo-molecular-simulation-software>`_.
Support FEASST with a `GitHub <https://github.com/usnistgov/feasst>`_ star or `manuscript <https://doi.org/10.1063/5.0224283>`_ citation!

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

How to compile
===============================================

FEASST requires C++14 and CMake, and is compiled with the following BASH commands:

.. code-block:: bash

    # [apt/yum/dnf/brew] install g++ cmake curl tar. On HPC, try "module avail/load"
    cd $HOME # replace this with your preference throughout
    curl -OL https://github.com/usnistgov/feasst/archive/refs/tags/v0.25.9.tar.gz # download
    tar -xf v0.25.9.tar.gz           # uncompress
    mkdir feasst-0.25.9/build; cd $_ # out-of-source build
    cmake ..                         # find prerequisites
    make install -j 4                # compile on 4 threads
    # Optional Python packages used in tutorials. Virtual environment recommended:
    # https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments
    pip install jupyter matplotlib pandas scipy ../pyfeasst # ../ ensures pyfeasst matches feasst
    # Because online documentation changes with verison, open local version documentation.
    xdg-open ../html/index.html      # on macOS, replace "xdg-open" with "open"

How to run a simulation
===============================================

Input a text file to the compiled executable.

.. code-block:: bash

    $HOME/feasst-0.25.9/build/bin/fst < $HOME/feasst-0.25.9/tutorial/example.txt

The following text input file is explained in detail in the first :doc:`tutorial <tutorial/tutorial>`.

.. literalinclude:: tutorial/example.txt


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

Troubleshooting install
================================================================================

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
