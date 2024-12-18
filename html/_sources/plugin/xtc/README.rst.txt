*********
XTC
*********

This plugin is used for reading and writing the XTC file format.

For writing compressed XTC trajectory files, you must first install the xdrfile library.

.. code-block:: bash

    wget ftp://ftp.gromacs.org/pub/contrib/xdrfile-1.1.tar.gz
    tar -xf xdrfile-1.1.tar.gz; cd xdrfile-1.1b
    mkdir build
    ./configure --enable-shared --prefix=$HOME/software/xdrfile-1.1b/build #enable-shared for SWIG
    make install
    export LD_LIBRARY_PATH="$HOME/software/xdrfile-1.1b/build/lib:$LD_LIBRARY_PATH"

Note, there is a newer XTC library but it leads to compiler warnings.

Associated CMake flag

.. code-block:: bash

   cmake -DXDRFILE_DIR=/path/to/xdrfile/build/ ..

.. toctree::
   :maxdepth: 1

   doc/toc

FEASST plugin dependencies
============================

* system
