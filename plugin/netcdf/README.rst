*********
NetCDF
*********

This is a plugin that enables an interface with NetCDF (https://www.unidata.ucar.edu/software/netcdf/).

The following are instructions for the manual installiation of zlib, hdf5 and netcdf.
You could also try to use your package manager.
Note the version numbers may change.

.. code-block:: bash

    # following https://support.scinet.utoronto.ca/education/staticpublic/course177content327.html
    mkdir $HOME/local
    mkdir $HOME/local/bin
    mkdir $HOME/local/include
    mkdir $HOME/local/lib
    # Add the next four lines to your .bashrc
    export PATH="$PATH:$HOME/local/bin"
    export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$HOME/local/lib"
    export LIBRARY_PATH="$LIBRARY_PATH:$HOME/local/lib"
    export CPATH="$CPATH:$HOME/local/include"
    cd ~/software
    wget http://zlib.net/zlib-1.2.13.tar.gz
    tar axvf zlib-1.2.13.tar.gz
    cd zlib-1.2.13
    ./configure --prefix=$HOME/local
    make
    make check
    make install
    cd ..
    wget https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.10.5.tar.gz
    tar axvf hdf5-1.10.5.tar.gz
    cd hdf5-1.10.5
    ./configure --prefix=$HOME/local
    make
    make check
    make install
    cd ..
    wget https://downloads.unidata.ucar.edu/netcdf-c/4.9.2/netcdf-c-4.9.2.tar.gz
    tar axvf netcdf-c-4.9.2.tar.gz
    cd netcdf-c-4.9.2/
    sudo apt install m4 curl
    ./configure --prefix=$HOME/local --enable-netcdf-4 --disable-libxml2 --disable_byterange
    make
    make check
    make install
    cd ..
    wget https://github.com/Unidata/netcdf-cxx4/archive/v4.2.1.tar.gz
    tar axvf v4.2.1.tar.gz
    cd netcdf-cxx4-4.2.1/
    ./configure --prefix=$HOME/local
    make
    make check
    make install
    cd ..

Now that the prerequisites are taken care of, FEASST must be installed again with some modification.
First, modify /path/to/feasst/CMakeLists.txt to include netcdf in set(FEASST_PLUGINS ...).
Also, include "-DUSE_NETCDF=ON -DNETCDF_DIR=$HOME/local" for the cmake command.

.. toctree::
   :glob:

   tutorial/tutorial*

FEASST plugin dependencies
============================

* system

API
===

.. toctree::
   :maxdepth: 1

   doc/toc
