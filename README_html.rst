*************************
README
*************************

The Free Energy and Advanced Sampling Simulation Toolkit (FEASST) is a free,
open-source, public domain software to conduct molecular- and particle-based
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

How to compile
===============================================

.. code-block:: bash

    #[apt/yum/dnf/brew] install g++ cmake git python3
    cd $HOME # replace this with your preference throughout
    git clone https://github.com/usnistgov/feasst.git
    mkdir feasst/build; cd $_
    cmake ..
    make install -j$(nproc)
    # optional python packages for feasst tutorials
    pip install jupyter matplotlib pandas scipy ../pyfeasst

FEASST requires a C++14 compiler, CMake and Python3, while git is optional.
The executables ``fst`` and ``rst`` should appear in ``$HOME/feasst/build/bin/``.
Text input files are run using ``fst < input.txt`` and simulations are restarted using ``rst checkpoint.txt``.
For pyfeasst, provide pip a path to the specific pyfeasst directory in feasst to ensure the versions match (e.g., do not leave out the "../").

Basic simulation example
===============================================

The following text input file is explained in detail in the first :doc:`tutorial <tutorial/tutorial>`.

.. literalinclude:: tutorial/example.txt

How to get started
===============================================

* Complete the :doc:`first <tutorial/tutorial>` and :doc:`second tutorial<tutorial/launch>`.

  * Copy/paste or use the URL to find the code (e.g., https://pages.nist.gov/feasst/tutorial/launch.html is ``$HOME/feasst/tutorial/launch.py``).
  * See ``python launch.py --help`` (e.g., adjust ``--feasst_install`` or ``--hours_terminate``).
* Find :doc:`tutorial/README` that are closest to what you would like to accomplish.
* Reproduce the expected result of those :doc:`tutorial/README`.
* To modify the tutorial to accomplish your goals, refer to the :doc:`../plugin/text_interface` documentation.
* Compare the energy of a :doc:`reference configuration<plugin/monte_carlo/tutorial/tutorial_0_ref_configs>` with a trusted source.

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
