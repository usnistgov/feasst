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

.. note::

   Website: https://feasst.hhatch.com

Features
================================================================================

* NVT/muVT LJ (with LRC) simulations with Metropolis or Wang-Landau acceptance criteria.
* Second virial coefficient using Mayer sampling for LJ with HS reference.

Compile
================================================================================

.. code-block:: bash

    mkdir /path/to/feasst/build
    cd /path/to/feasst/build
    cmake ..
    make -j12
    make test         # optional test

C++ Library Usage
================================================================================

.. code-block:: bash

    cd /path/to/feasst/build
    make install -j12
    cd /path/to/feasst/plugin/core/tutorial/
    mkdir build; cd build
    cmake ..
    make
    ./main

CMake defaults to install in the build directory.
But you can also specify the path as follows.

.. code-block:: bash

    cmake -DCMAKE_INSTALL_PREFIX=/path/to/install/dir ..

Later when you build your main executable you need to specify this path:

.. code-block:: bash

    cmake -DCMAKE_PREFIX_PATH=/path/to/install/dir ..

Python usage
================================================================================

.. code-block:: bash

    cd /path/to/feasst/build
    cmake -DUSE_SWIG=ON ..
    make _feasst -j12
    ../py/run.sh ../py/test.py

CMake attempts to find the python libraries.
But you may want to specify them manually as follows:

.. code-block:: bash

    cmake -DUSE_SWIG=ON -DSET_PYTHON_PATH=ON -DPYTHON_INCLUDE_PATH=/path/to/anaconda/include/python3.6m -DPYTHON_LIBRARIES=/path/to/anaconda/lib/libpython3.6m.so ..

.. include:: CONTACT.rst

.. include:: DISCLAIMER.rst

.. include:: LICENSE.rst
