*********
FFTW
*********

This is a plugin that enables an interface with the FFTW software (https://www.fftw.org/).

Installation

.. code-block:: bash

    cd $HOME/software/
    wget https://www.fftw.org/fftw-3.3.10.tar.gz
    tar -xf fftw-3.3.10.tar.gz
    #./configure --enable-shared --prefix=$HOME/software/fftw-3.3.10/build
    #./configure --enable-shared
    sudo make install -j24
    sudo apt-get install pkg-config # for cmake

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
