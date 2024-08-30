***************************************
MPI
***************************************

Use MPI with FEASST.

Setup prerequisites:

* sudo apt install openmpi-bin libopenmpi-dev
* sudo dnf install openmpi-devel; module load mpi/openmpi-x86_64
* sudo dnf install python3-devel python3.12-devel
* pip install mpi4py
* cmake -DUSE_MPI=ON ..
* recompile FEASST

.. toctree::
   :glob:

   tutorial/tutorial*

FEASST plugin dependencies
============================

* monte_carlo

API
===

.. toctree::
   :maxdepth: 1

   doc/toc
