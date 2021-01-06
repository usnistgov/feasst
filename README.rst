*************************
README
*************************

.. .. image:: https://anaconda.org/hhatch/feasst/badges/installer/conda.svg
..     :target: https://conda.anaconda.org/hhatch

.. image:: https://travis-ci.com/hhatch/feasst.svg?branch=main
    :target: https://travis-ci.com/hhatch/feasst

.. image:: https://codecov.io/gh/hhatch/feasst/branch/main/graph/badge.svg
    :target: https://codecov.io/gh/hhatch/feasst

.. image:: https://www.codefactor.io/repository/github/hhatch/feasst/badge/main
    :target: https://www.codefactor.io/repository/github/hhatch/feasst/overview/main
    :alt: CodeFactor

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

Simulation techniques

* Metropolis Monte Carlo
* Wang-Landau Monte Carlo
* Transition-matrix Monte Carlo
* Mayer-sampling Monte Carlo
* Configurational bias

Thermodynamic ensembles

* Microcanonical ensemble
* Canonical ensemble
* Grand canonical ensemble
* Expanded ensembles

Intermolecular interactions

* Hard spheres
* Lennard-Jones with Yukawa, LRC, cut and force shift
* Patchy particles
* Charged interactions with the Ewald summation
* Confinement
* Cell list and neighbor list

Modern software

* Interface with C++ or as a Python module
* OpenMP parallelization
* Checkpointing to save and restart simulations

Tutorial
================================================================================

`Tutorial <tutorial/tutorial.html>`_ describes an example Lennard-Jones Monte Carlo simulation.
After completing this basic tutorial, check out tutorials specific to each :doc:`/plugin/README`.
For example, see the tutorials of `System <plugin/system/README.html>`_, `MonteCarlo <plugin/monte_carlo/README.html>`_ and `FlatHistogram <plugin/flat_histogram/README.html>`_.

The search box for the html documentation was disabled for security reasons.
But the html is generated entirely from the downloaded code.
Thus, the grep command is a great option to search for more information on keywords.
For example, if you would like more information on `RandomMT19937 <plugin/math/doc/RandomMT19937.html>`_ but are not sure where to find it, you could search headers files

.. code-block:: bash

   grep -r --include=*.h RandomMT19937

And find that the class is part of the `Math <plugin/math/README.html>`_ plugin.

Build from source code
=======================

To begin, [apt/yum/brew] install g++ cmake git.

Python install
----------------

* SWIG is required. Version 3.0.12 is recommended if your current SWIG version does not work properly.

* Anaconda with Python 3 is recommended.

* CMake attempts to find the python libraries during compilation.
  But you may want to specify them manually.

.. code-block:: bash

    git clone https://github.com/usnistgov/feasst.git
    mkdir feasst/build
    cd feasst/build

    conda env create -f ../py/feasst.yml
    conda activate feasst
    # alternatively, without env create, include some tutorial dependencies
    # pip install jupyter matplotlib pandas scipy # for tutorials

    cmake -DUSE_SWIG=ON ..
    # alternatively, for manually setting the python path
    # cmake -DUSE_SWIG=ON -DSET_PYTHON_PATH=ON -DPYTHON_INCLUDE_DIR=/path/to/anaconda/include/python3.7m -DPYTHON_LIBRARY=/path/to/anaconda/lib/libpython3.7m.so ..

    make -j4
    make install -j4
    python ../py/test.py # optional test

C++ install
----------------

First, install the C++ library.

.. code-block:: bash

    git clone https://github.com/usnistgov/feasst.git
    mkdir feasst/build
    cd feasst/build
    cmake ..          # optionally, include -DUSE_GTEST=ON for gtest
    make install -j4
    make test         # optional test

Then, compile the specific simulation you wish to run (e.g., tutorial).

.. code-block:: bash

    cd /path/to/feasst/tutorial/
    mkdir build; cd build
    cmake ..
    make
    ./tutorial

CMake defaults to install in the build directory.
But you can also specify the path as follows.

.. code-block:: bash

    cmake -DCMAKE_INSTALL_PREFIX=/path/to/install/dir ..

Later, when you build your tutorial executable, if your build directory is not `~/feasst/build`, then specify the path to the build directory as follows:

.. code-block:: bash

    cmake -DCMAKE_PREFIX_PATH=/path/to/install/dir ..

Sometimes running the executable results in an error that a plugin cannot be found.
In that case `LD_LIBRARY_PATH` needs to be set as follows:

.. code-block:: bash

    export LD_LIBRARY_PATH="/path/to/feasst/build/:$LD_LIBRARY_PATH

Troubleshooting install
------------------------

Please :doc:`/CONTACT` us if you run into an issue not listed below.

macOS Mojave
~~~~~~~~~~~~~~

* SWIG (from Homebrew) is likely version 4, which sometimes causes a SEGFAULT when trying to run feasst.
  Try SWIG version 3 instead.

* Sometimes CMake has trouble finding anaconda, and if you use SET_PYTHON_PATH described above, you may need to look out for the .dylib instead of .so

CentOS 7
~~~~~~~~~

CMake and SWIG versions are usually too old.
Try the command cmake3 instead of cmake.
Otherwise, install SWIG 3.

Windows 10
~~~~~~~~~~~

* Install Windows subsystem for Linux (Ubuntu 16)
* See Ubuntu 16

Ubuntu 16
~~~~~~~~~~

* Update to CMake 3 (https://cmake.org/download/)
* sudo apt install swig

Ubuntu 18 or 20
~~~~~~~~~~~~~~~~

* We are not aware of any issues with an Ubuntu 18 or 20 install.

Cray (NERSC CORI)
~~~~~~~~~~~~~~~~~~

* OpenMP functions apparently do not work unless the cray programming environment is disabled.

Build from Docker
===================

The Docker is no longer being updated and contains a very old version of the code.
Installation via `Docker <docker.io>`_ is an alternative if you are experiencing installation issues.
Unfortunately there is a performance penalty using this method.

First, install docker

 * Ubuntu: `apt install docker`
 * macOS: download and install Docker Desktop.
 * CentOS: `yum install docker`

On a mac, run Docker.app, login, startup and skip this step.
Otherwise, start docker and create a `docker group <https://docs.docker.com/install/linux/linux-postinstall/>`_ to avoid root.

.. code-block:: bash

    sudo service docker start
    sudo groupadd docker
    sudo usermod -aG docker $USER
    newgrp docker

.. code-block:: bash

    git clone https://github.com/usnistgov/feasst.git
    docker build -t hhatch/feasst:firsttry feasst/

Alternatively you could pull feasst from docker hub (version may be old):

.. code-block:: bash

    docker pull hhatch/feasst:firsttry

Then you can interactively run feasst via

.. code-block:: bash

    docker run -v ~/scripts:/mnt/scripts -it hhatch/feasst:firsttry
    cd /mnt/scripts

Or with `~/notebooks`

.. code-block:: bash

    docker run -v ~/notebooks:/mnt/notebooks -i -t -p 8888:8888 hhatch/feasst:firsttry /bin/bash -c "jupyter notebook --notebook-dir=/mnt/notebooks --ip='*' --port=8888 --no-browser --allow-root"

.. include:: CONTACT.rst

.. include:: DISCLAIMER.rst

.. include:: LICENSE.rst
