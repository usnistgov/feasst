**********************************************************************
FEASST - Free Energy and Advanced Sampling Simulation Toolkit
**********************************************************************
Author: Harold Wickes Hatch, hhatch.com, harold@hhatch.com

FEASST is a free, open-source program to conduct molecular and 
particle-based simulations with flat-histogram Monte Carlo methods.

Features
- Wang-Landau, Transition-Matrix and/or Metropolis Monte Carlo
- Canonical, grand canonical and expanded ensembles
- Interface as a Python module or C++ class
- Easy use of MPI and checkpointing (stop and restart simulations)
- Many advanced sampling Monte Carlo moves
- Supported energy functions include anisotropic particles and the
    Ewald summation
- Robust unit testing and thoroughly documents case studies

Advanced Monte Carlo moves
- Configurational bias insertions, deletions and regrowth with 
    multiple first bead insertion
- Aggregation volume bias (AVB) insertions, deletions and the AVB2 
    and AVB3 algorithms
- Geometric cluster algorithm
- Rigid cluster moves
- Parallel configuration swaps

Intermolecular interactions
- Charged interactions with the Ewald summation
- Lennard Jones with Yukawa, LRC, force shift, gaussians
- Patchy particles
- Hard spheres/square well
- Tabular potentials for anisotropic particles

Documentation may be in html and latex directories (requires Doxygen)

**********************************************************************
Requirements
**********************************************************************

FEASST is designed for a LINUX or MAC OS X platform.

make 3.81
g++ 4.7 (support for c++0x standard)
xdrfile 1.1b (included)

Optional tools:
gtest-1.7.0 (for unittests)
valgrind (optional for development)
doxygen 1.6.1 (optional for documentation)
swig 1.3.40 (optional for python interface)
anaconda 1.9.1 (optional, recommended for python compatibility)
openmpi 1.4.5 (optional)

**********************************************************************
Installation
**********************************************************************

For basic installation:
  cd src
  make cnotest
  #NOTE: if you are not using xtc files, remove the xdr flags in LIBS
  # and CXXFLAGS in src/Makefile

FEASST needs to know where its install directory is located.
  This is done by adding the following command in your ~/.bash_profile
  export FEASST_INSTALL_DIR_="/put/your/path/to/feasst/here/"
  This is important for run.sh and compile.sh

For XTC files:
  ftp://ftp.gromacs.org/pub/contrib/xdrfile-1.1.tar.gz
  tar -xf xdrfile-1.1.tar.gz; cd xdrfile-1-1b
  ./configure --enable-shared --prefix=$HOME/ #enable-shared for swig
  make install
  #NOTE: You may install xdrfile in a different install location
  # via prefix, which requires you to include the path in CXXFLAGS
  # of src/Makefile

For unittests
  download gtest: http://code.google.com/p/googletest/downloads/detail?name=gtest-1.7.0.zip
  unzip and correct path in GTEST_DIR in Makefile. No need to compile.
  to compile the FEASST unittests, use the command 'make c' in src/

Download and install openmpi 1.4.5 with Intel compilers
  tar -xf openmpi*gz; cd openmpi*; mkdir build; cd build
  ../configure --prefix=`pwd`/.. CC=icc CXX=icpc 
  make
  make install

FFTW3.3.4 (optional)
  download fftw-3.3.4, uncompress, move to main directory
  ./configure --prefix=/path/to/install/dir --enable-shared --with-pic
  make
  make install

VMD (optional)
  download vmd
  tar -xf vmd-1.9.2.bin.LINUXAMD64-RHEL5.opengl.tar.gz
  cd vmd-1.9.2
  edit the configure file to change install location
  ./configure LINUXAMD64
  make install -j 8
  add VMD to your path
   e.g. , in your .bash_rc file, 
   PATH=$PATH:$HOME/software/vmd-1.9.2/bin/

SWIG 2.0.12 (optional)
  cd swig-2.0.12; ./configure --prefix=/path/to/dir; make

HDF5-1.8.18 (optional)
  sudo ./configure --prefix=/usr/local/hdf5 --enable-cxx
  make; make check; make install; make check-install
  ##add library to LD_LIBRARY_PATH

DO NOT USE tools/compile.sh FOR INITIAL INSTALLATION

**********************************************************************
Usage and Test Cases
**********************************************************************

FEASST may be accessed through a C++ or Python interface.

For the C++ interface, the compile.sh script in the tools directory
is used to compile the C++ interface file to generate an executable
for simulation.

For the Python interface, the path to FEASST is on line 5

**********************************************************************
Test Case - Lennard Jones
**********************************************************************

The first test case is to simulate the LJ potential and reproduce the
results found in the NIST SRSW for LJ, EOS-TMMC: 
http://www.nist.gov/mml/csd/informatics_research/srsw.cfm
http://mmlapps.nist.gov/srs/LJ_PURE/eostmmc.htm

Move to the test directory
cd ../test/lj/srsw/eostmmc/1.5

The file muvttmmclj.cc is the C++ interface for FEASST.
The file run.sh is an example script to compile and run muvttmmclj
The file data.lj is a LAMMPS style description of an LJ particle
The files lj.msdb.* contain the results that you may compare with

Run the simulation: ./run.sh, or ./muvttmmclj.py, and it will produce:
log - a log file printing information by step
colMat - a file with the macrostate probability distribution (lnpi),
  potential energy and collection matrix
movie* - movie files viewable with VMD

To compare the results with the NIST SRSW, compare the following:
  colMat columns 1:2 with lj.msdb.t150.*.p_macro.dat
  colMat columns 1:3 with lj.msdb.t150.*.energy.dat

To run the test on multiple processors using OMP:
  cd omp
  ./run.sh OR ./muvttmmclj.py

**********************************************************************
Example of a 1D Tabular Potential
**********************************************************************

Instead of implementing your own pair potential in the code itself, 
you may simply make a text file with your potential. An example may be
found in the following test directory:

test/binary/tee/table

In this example, a binary LJ-lambda potential is simulated. In tee.cc
the potential is implemented with PairLJSlow, printed, and then used
to initialize PairTabular. It outputs the tables as files "tabi*j*"
which have headers like the following:

  # tableType PairLJSlow
  # tabDims 1
  # dimorder0 0
  # dimn0 5001
  # dimmin0 0.94089999999999996
  # dimmax0 1.1664000000000001
  9542.2200121376991
  9483.2587236908766
  9424.6627162728782

Note that dimn0 is the number of table elements. Distances are shown
as a function of the variable s=r^2, such that dimmax0=rCut^2=1.08^2
and dimmin0=rCutInner^2=sigFac^2=0.97^2. For tabular potentials, 
r < rCutInner has infinite potential energy.

An example of utilzing the table potential (without generating) is
provided in the following file:

test/binary/tee/table/tee_nogen.cc

This file and the tabi*j* files may be used as templates to create
your own pair potentials.

**********************************************************************
Analysis of Configurations for WLTMMC simulations
**********************************************************************

For analyzing configurations to post-process simulation data (e.g., 
log and movie files), see the example code:

tools/xyz2bin.cc

This code both splits the xyz files and shows example of doing some 
analysis within the c code. It uses the checkpointing to pick up run 
variables. For example it knows how often you printed movies versus 
logs. It also picks up on the number of processors and can average 
over all of those processors for analysis.

It outputs files:
  analysis* for average x position of first molecule
  moviep[proc]b[bin].xyz

where proc is the processor number and bin is the order parameter
index as descripted by the acceptance criteria.

**********************************************************************
Initializing a Simulation from an XYZ File
**********************************************************************

The following code reads an xyz file format to input an initial 
configuration.

  std::ifstream inFile("nameOfFile");
  p.readxyz(inFile);  // "p" is the name of your pair object

If there are zero particles in the space class, then it automatically 
attempts to add the molecules based on the first molecule type 
described by the s.addMolInit function.

In cases with multiple components, this is not sufficient. So you will
want to make sure you initialize the appropriate number of molecules 
in the appropriate order. For example, if your xyz file lists A,B,A,B.
or A,A,A...,B,B,B... then you need to add these in the right order.

For example, something like the following:

for (int iMolA = 0; iMolA < nMolA; ++iMolA) {
  s.addMol("/name/of/data/file");
}

Then the same for B, assuming your xyz has all A listed, followed by 
all B.

After all of the s.addMol commands are performed and the xyz file is 
read, you will need to update the pair class as follows

  p.addPart();
  p.initEnergy();

A simple test that the xyz file was read correctly is to print it and
compare:

  p.printxyz("filename",1);

**********************************************************************
Example of Restarting a Simulation
**********************************************************************

Checkpoint files may be written periodically during a simulation, and
these may be used to restart a simulation. For example, see the file

test/binary/tee/table/tee_rst.cc

In this file, the checkpoint file is simply read and restarted in
two lines for single processor simulations:

  // read checkpoint files
  WLTMMC mc("tmp/rst");

  // run simulation
  mc.runNumTrials(npr);

Note that, if you are attempting to restart a simulation that was
terminated abruptly, it is possible that the checkpoint files were in
the process of writing during the termination. In this case, the files
themselves could be missing important details. If this is the case,
your simulation will likely crash upon restart or output potential
energy which is quite different from the previous value. To remedy
this situation, the checkpoint files keep a 'backup' file which ends
in ".bak", which you may use instead. If you wish to use the backup
files, then all files ending in '.bak' should replace the same files
without the '.bak' ending. Before replacing files, it is recommended
to first backup the entire tmp folder.

Note that multiprocessor simulations may take additional care to
restart correctly. If you wish to restart just one processor, you may
simply use one of the files with an appended p# (e.g. "tmp/rstp0").

If you wish to restart simulations that are independent, then an 
example may be found in the following file

/home/hwh/feasst/test/binary/tee/shortrange/tee_rst.cc

In this file, the two lines are as follows:

  // read restart file
  WLTMMC mc("tmp/rst");

  // run sweeps
  mc.runNumSweepsRestart(100, "tmp/rst");

Restarting simulations that are coupled (e.g. by configuration swaps)
may require more initialization that is not currently described in 
this documentation.

************************************************************
For more advanced documentation, see DEV.txt for developers
************************************************************
