*********
FFTW
*********

This is a plugin that enables an interface with the FFTW software (https://www.fftw.org/).

Installation of prerequisites

.. code-block:: bash

    cd $HOME/software/
    wget https://www.fftw.org/fftw-3.3.10.tar.gz
    tar -xf fftw-3.3.10.tar.gz
    cd fftw-3.3.10
    ./configure --enable-shared
    sudo make install -j24
    sudo apt-get install pkg-config # for cmake to find fftw

Now that the prerequisites are taken care of, FEASST must be installed again with some modification.
First, modify /path/to/feasst/CMakeLists.txt to include fftw in set(FEASST_PLUGINS ...).
Also, include -DUSE_FFTW=ON for the cmake command.

.. toctree::
   :glob:

   tutorial/tutorial*

FEASST plugin dependencies
============================

* system

API
===

.. toctree::

   doc/toc
