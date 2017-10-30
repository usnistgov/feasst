**********************************
Frequently asked questions (FAQ)
**********************************

How do I get started?
#######################

Try running some test cases and find one that is closest to what you want to do.
Reproduce a published result and carefully modify it as you see fit.

Simulation issues
###################################################

Calculate average properties
=========================================

* The log file from :cpp:func:`MC::initLog()` may periodically output the instantaneous value of the property of interest.
* Obtain the property of interest from the :doc:`/api` periodically, as shown in the python script of :doc:`/testcase/1_lj_2_nvt-mc_README`.
* Implement a custom Analysis class as demonstrated for average energy in the C++ script of :doc:`/testcase/1_lj_2_nvt-mc_README`

Calculate grand canonical ensemble average properties
=================================================================

* See the example in :doc:`/testcase/1_lj_3_eostmmc_README` for coexistence properties.

Load coordinate files or place particles in specific locations
================================================================

* See the example in :doc:`/testcase/1_lj_1_ref-config_README`.

Make FEASST do exactly what you want it to do
===========================================================

You can create your own :cpp:class:`Pair`, :cpp:class:`Trial`, :cpp:class:`Criteria`, :cpp:class:`Analyze`, and :cpp:class:`Random` classes as described in the documentation for the Base classes.
Start by finding an existing derived class that is the most similar and copying it.
You can also define new classes within the script themselves.
For example, see :doc:`/testcase/1_lj_2_nvt-mc_README` for a custom :cpp:class:`Analyze`, and :doc:`/testcase/2_jagla_1_ref-config_README` for a custom :cpp:class:`Pair`.

Analysis of configurations for WL-TMMC simulations
==================================================================================

For analyzing configurations to post-process simulation data (e.g.,
log and movie files), see the example code `<tools/xyz2bin.cc>`_

This code both splits the xyz files for a given order parameter and shows example of doing some
analysis within the c code. It uses the checkpointing to pick up run
variables. For example it knows how often you printed movies versus
logs. It also picks up on the number of processors and can average
over all of those processors for analysis.

It outputs files::

    `analysis*` for average x position of first molecule
    moviep[proc]b[bin].xyz

where proc is the processor number and bin is the order parameter
index as described by the acceptance criteria.

Initializing a simulation from an XYZ file
============================================

The following code reads an xyz file format to input an initial
configuration.

.. code-block:: C++

    std::ifstream inFile("nameOfFile");
    pair.readxyz(inFile);

If there are zero particles in the space class, then it automatically
attempts to add the molecules based on the first molecule type
described by the s.addMolInit function.

In cases with multiple components, this is not sufficient. So you will
want to make sure you initialize the appropriate number of molecules
in the appropriate order. For example, if your xyz file lists A,B,A,B.
or A,A,A...,B,B,B... then you need to add these in the right order.

For example, something like the following:

.. code-block:: C++

    for (int iMolA = 0; iMolA < nMolA; ++iMolA) {
      pair.addMol("/name/of/data/file");
    }

Then the same for B, assuming your xyz has all A listed, followed by
all B.

After all of the s.addMol commands are performed and the xyz file is
read, you will need to update the pair class as follows:

.. code-block:: C++

    p.initEnergy();

A simple test that the xyz file was read correctly is to print it and
compare:

.. code-block:: C++

    p.printxyz("filename",1);

Restarting a simulation
=========================

Checkpoint files may be written periodically during a simulation, and
these may be used to restart a simulation. For example, see `<test/binary/tee/table/tee_rst.cc>`_

In this file, the checkpoint file is simply read and restarted in
two lines for single processor simulations:

.. code-block:: c++

    // read checkpoint files
    feasst::WLTMMC mc("tmp/rst");

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
example may be found in `<tools/rstMultiproc.cc>`_.

In this file, the two lines are as follows:

.. code-block:: c++

    // read restart file
    feasst::WLTMMC mc("tmp/rst");

    // run sweeps
    mc.runNumSweepsRestart(100, "tmp/rst");

Restarting simulations that are coupled (e.g. by configuration swaps)
may require more initialization that is not currently described in
this documentation.

Isotropic tabular potential
=============================

Instead of implementing your own pair potential in the code, you may simply make a text file with your potential.
An example may be found in the following test directory: `<test/binary/tee/table>`_

In this example, a binary LJ-lambda potential is simulated. In tee.cc
the potential is implemented with PairLJMulti, printed, and then used
to initialize PairTabular. It outputs the tables as files ``tabi*j*``
which have headers like the following::

    # tableType PairLJMulti
    # tabDims 1
    # dimorder0 0
    # dimn0 5001
    # dimmin0 0.94089999999999996
    # dimmax0 1.1664000000000001
    9542.2200121376991
    9483.2587236908766
    9424.6627162728782

Note that ``dimn0`` is the number of table elements.
Distances are shown as a function of the variable s=:math:`r^2`, such that dimmax0 = :math:`rCut^2 = 1.08^2` and dimmin0 = :math:`rCutInner^2 = sigFac^2 = 0.97^2`.
For tabular potentials, r < rCutInner has infinite potential energy.

An example of utilizing the table potential (without generating) is provided in
`<test/binary/tee/table/tee_nogen.cc>`_

This file and the ``tabi*j*`` files may be used as templates to create
your own pair potentials.

Getting to high density
==========================

This is a common issue with a few approaches:

1. You can use the parallel version with multiple windows. The ``MC::initWindows(nExp=`` variable sets the spacing based on some exponential distribution,e.g., nExp=3. makes windows even bigger o the low density side vs nExp=2.

http://feasst.hhatch.com/WLTMMC.html#_CPPv2N6WLTMMC11initWindowsEKdKi

2. Break into multiple simulations. For example, I launched two jobs simultaneously, one with N=0 to 336, and another N=336 to N=550. Of course the low density ones finished 10x faster but then the high density windows were smaller

3. If you're trying to do a high throughput approach where each model and condition can be quite different but you don't want to have to hand pick conditions for each one, you can do a ridiculous number of windows (N=16/32 perhaps?) and terminate after some run time and only use the windows that managed to converge and throw out the rest. Still you would want to choose carefully the nMolMax to not have too much wasted processor time.

4. You can use the more MD-style approach of initializing a big box and squeezing it.

.. code-block:: c++

   c.pressureset(1000.);
   transformTrial(&mc, "lxmod", 0.001);
   transformTrial(&mc, "lymod", 0.001);
   transformTrial(&mc, "lzmod", 0.001);
   // transformTrial(&mc, "vol", 0.001);  // this one does an isotropic volume move instead of independent. You don't need both

To remove the box move after equilibration, I prefer to "scope" the mc class with the box move and just make another clean mc class, or do an shallow copy (cloneShallow) before calling the transform trial and then use that clone later for production.
Another option is mc.removeTrial(trial#);

Compilation issues
###################################################

ModuleNotFoundError: "No module named 'feasst'"
================================================

* Verify that "/path/to/feasst/build/_feasst.so" was built without error, as described in :doc:`/readme`.
* Use the command "`<tools/run.sh>`_ script.py" instead of "``python script.py``"
* Alternatively, point to build/_feasst.so in your python script as follows

.. code-block:: py

   import sys
   sys.path.append(/path/to/feasst/build/)

error: ‘[function/class]’ is not a member of ‘feasst’
================================================================================================

In C++, this error often occurs when one does not include a necessary header file.
For example, to use `feasst::swapTrial()`, one must have `#include "trial_swap.h"`.

Compilation error gives "can not be used when making a shared object; recompile with -fPIC"
================================================================================================

One of your external libraries (e.g., fftw or xdrfile) needs the flag "--enable-shared" during configuration.
Or you can edit CMakeLists.txt to add "-fPIC" as follows:

SWIG_LINK_LIBRARIES(feasst ${PYTHON_LIBRARIES} ${EXTRA_LIBS} -fPIC)   # HWH: add -fPIC
#SWIG_LINK_LIBRARIES(feasst ${PYTHON_LIBRARIES} ${EXTRA_LIBS})        # HWH: old version


