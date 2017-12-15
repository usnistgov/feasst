*************************
README
*************************

The Free Energy and Advanced Sampling Simulation Toolkit (FEASST) is a free,
open-source, modular program to conduct molecular and particle-based
simulations with flat-histogram Monte Carlo and molecular dynamics methods.

.. note::

   DOI: https://doi.org/10.18434/M3S095

   Website: https://pages.nist.gov/feasst/

   Code repository: https://github.com/usnistgov/feasst

Features
########

.. image:: sphinx/feasst.png
   :target: https://pages.nist.gov/feasst
   :align: right

Features with the (dev) label are under development and may not be available in your release.
Please `Contact`_ us if interested.

Simulation techniques

* Wang-Landau Monte Carlo
* Transition-matrix Monte Carlo
* Metropolis Monte Carlo
* Mayer sampling Monte Carlo (dev)
* Molecular dynamics (dev)

Thermodynamic ensembles

* Microcanonical ensemble
* Canonical ensemble
* (Semi-)grand canonical ensemble
* Expanded ensembles in temperature, shape, etc. (dev)

Advanced Monte Carlo moves

* Parallel configuration swaps
* Constant pressure and floppy box
* Particle identify and position swaps
* Configurational bias insertions, deletions and regrowth with multiple first
  bead insertion (dev)
* Aggregation volume bias (AVB) insertions, deletions and the AVB2 and AVB3
  algorithms (dev)
* Geometric cluster algorithm (dev)
* Rigid cluster moves (dev)

Intermolecular interactions

* Hard spheres, soft spheres and square wells
* Charged interactions with the Ewald summation
* Lennard-Jones with Yukawa, LRC, force shift, or Gaussian
* Superquadrics and supertoroids (dev)
* Patchy particles (dev)
* Cylindrical and slit pore confinement (dev)

Modern software

* Interface with C++ or as a Python module
* OMP parallelization
* Checkpointing to save and restart simulations
* Robust unit testing

Installation
#############

FEASST is designed for a LINUX or MAC platform with the following minimum version software.

* make >= 3.81
* CMake >= 2.8.12.2
* compiler with c++0x support (e.g., g++ >= 4.7)

.. code-block:: bash

    git clone https://github.com/usnistgov/feasst.git
    cd feasst
    mkdir build
    cd build
    cmake ..
    make

Usage: C++ interface
#######################

The following may be found in the `<testcase/1_lj/0_example>`_ directory.

In C++, a simple NVT Lennard-Jones (LJ) simulation is performed as follows:

.. literalinclude:: testcase/1_lj/0_example/test.cc
   :language: c++

This C++ code is compiled and run in bash as follows:

.. code-block:: bash

    $HOME/feasst/tools/run.sh test.cc

Usage: Python interface
#########################

Requirements

* SWIG >= 1.3.40
* anaconda >= 1.9.1 (python >= 2.7)

To install the python interface, use the following CMake command in place of "cmake ..":

.. code-block:: bash

    cmake -DUSE_SWIG=ON -DPYTHON_INCLUDE_PATH=/path/to/anaconda/include/python3.6m -DPYTHON_LIBRARIES=/path/to/anaconda/lib/libpython3.6m.so .
    make _feasst -j

Note that the ``PYTHON_INCLUDE_PATH`` and ``PYTHON_LIBRARIES`` depends on your python installation.

The following may be found in the `<testcase/1_lj/0_example>`_ directory.
In python, a simple NVT Lennard-Jones (LJ) simulation is performed as follows:

.. literalinclude:: testcase/1_lj/0_example/test.py
   :language: py

This simulation is run in bash as follows:

.. code-block:: bash

    $HOME/feasst/tools/run.sh test.py

Optional external libraries
#######################################

* xdrfile 1.1b (compressed xtc trajectories)
* gtest >= 1.7.0 (C++ unittests)
* valgrind (C++ memory testing for development)
* doxygen >= 1.6.1 (C++ documentation)
* openmpi >= 1.4.5 (parallel computation)

To control the install, you can edit ``CMakeLists.txt`` in ``build`` as follows
before running the ``cmake ..`` command.

To use the XDRFILE library for xtc files:

.. code-block:: cmake

    option(USE_XDRFILE "Use xdrfile library" ON)

Or

.. code-block:: bash

    cmake -DUSE_XDRFILE=ON .

To give CMake the path to your xdrfile library:

.. code-block:: cmake

    set(XDRFILE_DIR "/path/to/xdrfile")

Or

.. code-block:: bash

    cmake -DXDRFILE_DIR=/path/to/xdrfile .

If you are changing the default build options in ``CMakeLists.txt``,
make sure to start compilation with a fresh ``build`` directory before CMake is
invoked (e.g., completely remove the build directory and start over, after
saving any relevant changes to ``CMakeLists.txt``).

Here is how to set up external libraries you may want to use with FEASST.
To begin, some libraries require installation.

XTC 1.1b
********

For writing compressed XTC trajectory files.

.. code-block:: bash

    ftp://ftp.gromacs.org/pub/contrib/xdrfile-1.1.tar.gz
    tar -xf xdrfile-1.1.tar.gz; cd xdrfile-1-1b
    mkdir build
    ./configure --enable-shared --prefix=$HOME/software/xdrfile-1.1b/build #enable-shared for SWIG
    make install
    export LD_LIBRARY_PATH="$HOME/software/xdrfile-1.1b/build/lib:$LD_LIBRARY_PATH"

Associated CMake flag

.. code-block:: bash

   cmake -DUSE_XDRFILE=On -DXDRFILE_DIR=/path/to/xdrfile .

Google Test 1.7.0
*****************

For testing the C++ code: CMake automatically clones and compiles the repository.

Associated CMake flag

.. code-block:: bash

   cmake -DUSE_GTEST=On .

.. OpenMPI with Intel compilers
   ****************************
   .. code-block:: bash
       tar -xf openmpi*gz; cd openmpi*; mkdir build; cd build
       ../configure --prefix=`pwd`/.. CC=icc CXX=icpc $intel compilers
       make
       make install

OpenMP
******

CMake automatically searches for OpenMP support from the compiler.

FFTW 3.3.4
**********

This library is used for computing the scattering of anisotropic shapes.

.. code-block:: bash

    # download fftw-3.3.4, uncompress, move to main directory
    ./configure --prefix=/path/to/install/dir --enable-shared --with-pic
    make
    make install

Associated CMake flag

.. code-block:: bash

   cmake -DUSE_FFTW=On -DFFTW_DIR=/path/to/fftw .

VMD 1.9.2
*********

VMD is great for visualizing and analyzing trajectories.

.. code-block:: bash

    # download vmd
    tar -xf vmd-1.9.2.bin.LINUXAMD64-RHEL5.opengl.tar.gz
    cd vmd-1.9.2
    # edit the configure file to change install location
    ./configure LINUXAMD64
    cd src
    make install -j 8
    # add VMD to your path
    export PATH=$PATH:/path/to/install/dir/vmd-1.9.2/bin/
    # I've noticed on centos6 or rocks6, export LIBGL_ALWAYS_INDIRECT=yes

SWIG 2.0.12
************

Required for python installation.

.. code-block:: bash

    cd swig-2.0.12; ./configure --prefix=/path/to/install/dir; make; make install

Associated CMake flag

.. code-block:: bash

   cmake -DUSE_SWIG=On .

CMake 2.8.12.2
**************

Download from https://cmake.org/files/v2.8/ ::

    tar -xf cmake-2.8.12-rc2-Linux-i386.tar.gz

HDF5 1.8.18
***********

.. code-block:: bash

    sudo ./configure --prefix=/usr/local/hdf5 --enable-cxx
    make; make check; make install; make check-install

Associated CMake flag

.. code-block:: bash

   cmake -DUSE_HDF5=On -DHDF5_USER_DIR=/path/to/hdf5 .

GSL 2.3
*******

For spline interpolation.

.. code-block:: bash

    ./configure --prefix=/path/to/install/dir; make; make install

Associated CMake flag

.. code-block:: bash

   cmake -DUSE_GSL=On -DGSL_USER_DIR=/path/to/gsl .

LCOV 1.13-1
***********

Required for html output of CMake command ``make coverage``
For graphical front-end of gcov, http://ltp.sourceforge.net/coverage/lcov.php ::

    rpm -i lcov-1.13-1.noarch.rpm

Associated CMake flag

.. code-block:: bash

   cmake -DUSE_GCOV=On .

Contact
#######

Project lead: Harold Wickes Hatch

https://www.nist.gov/people/harold-hatch

harold.hatch@nist.gov

For list of contributors, see `<CONTRIBUTORS.rst>`_

.. include:: CITATION.rst

.. include:: DISCLAIMER.rst

.. include:: LICENSE.rst
