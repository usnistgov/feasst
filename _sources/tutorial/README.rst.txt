***************************
Tutorial
***************************

.. While this rst file is human-readable, it is optimized for generating the HTML tutorial pages.

Installation and Interface
===========================

.. sidebar:: Python or C++?

   While the majority of FEASST users prefer the Python interface, FEASST is written almost entirely in C++.
   Thus, both interfaces will be supported in the long term.

The installation is described in the :doc:`/README`.
Although this tutorial will focus upon the python interface, FEASST may also be utilized as a C++ library.

.. note::

   Although this tutorial focuses upon the Python interface, a comparison of `tutorial/lj_brief.py` with `tutorial/lj_brief.cpp` shows that the C++ and Python interfaces are nearly identical.
   Thus Python tutorials are sufficient for learning the C++ library.

.. _lj_brief_tutorial:

Canonical ensemble Lennard-Jones Monte Carlo
=======================================================

The following simulation (tutorial/lj_brief.py) demonstrates the basics of FEASST.

.. literalinclude:: lj_brief.py
   :language: py
   :linenos:

.. sidebar:: :cpp:class:`MonteCarlo <feasst::MonteCarlo>` class diagram

   Top to bottom and left-to-right order is important.
   For example, acceptance the :cpp:class:`Criteria <feasst::Criteria>` should be defined before the :cpp:class:`Trials <feasst::Trial>`.
   Solid arrows show composition while dashed arrows show inheritance.

:cpp:class:`MonteCarlo <feasst::MonteCarlo>` contains :cpp:class:`System <feasst::System>`, acceptance :cpp:class:`Criteria <feasst::Criteria>`, :cpp:class:`Trials <feasst::Trial>` and, finally, the usually infrequent :cpp:class:`Analyze <feasst::Analyze>` and :cpp:class:`Modify <feasst::Modify>` between :cpp:class:`Trials <feasst::Trial>`.
:cpp:class:`System <feasst::System>` contains :cpp:class:`Configurations <feasst::Configuration>` and :cpp:class:`Potentials <feasst::Potential>`.
:cpp:class:`Configurations <feasst::Configuration>` contain :cpp:class:`Particles <feasst::Particle>` and the :cpp:class:`Domain <feasst::Doman>`.
The interaction :cpp:class:`Potentials <feasst::Potential>` contain :cpp:class:`Models <feasst::Model>` and :cpp:class:`VisitModels <feasst::VisitModel>`.
The acceptance :cpp:class:`Criteria <feasst::Criteria>` may be either :cpp:class:`Metropolis <feasst::CriteriaMetropolis>` or :cpp:class:`flat histogram <feasst::CriteriaFlatHistogram>`.
This dependency is summarized in the following Figure and is useful to remember when initializing any :cpp:class:`MonteCarlo <feasst::MonteCarlo>` simulation.

.. figure:: ../plugin/monte_carlo/doc/monte_carlo.svg

Thus, a first step in initializing :cpp:class:`MonteCarlo <feasst::MonteCarlo>` is to add a :cpp:class:`Configuration <feasst::Configuration>` (or a :cpp:class:`System <feasst::System>`).

.. literalinclude:: lj_brief.py
   :language: py
   :lines: 4-6
   :lineno-match:

.. sidebar:: :cpp:class:`Argument <feasst::Arguments>` types

   Although some type checking is built into :cpp:class:`Arguments <feasst::Arguments>`, care must be taken to input strings which follow the documentation of the class, as found in the :doc:`/plugin/README`.

In this example, a simple cubic periodic box of length of 8 is defined with a single type of particle as described in `forcefield/data.lj`.
See :doc:`/forcefield/README` for more information about the format of the data file, which is a LAMMPS-inspired file with some major differences.
:cpp:class:`Arguments <feasst::Arguments>` are input by a dictionary of strings.
The SWIG C++ Python interface converts this Python dictionary of strings into a C++ std::map of strings via `feasst.args(...)`.

Initialization of the :cpp:class:`Potentials <feasst::Potential>` proceeds as follows:

.. literalinclude:: lj_brief.py
   :language: py
   :lines: 7-8
   :lineno-match:

.. sidebar:: A FEASST convention

  :cpp:func:`MakeModelLJ <feasst::MakeModelLJ>` creates a new :cpp:class:`ModelLJ <feasst::ModelLJ>` object and returns a pointer to that object.
  This serves two purposes involving C++11 smart pointers and brace enclosed initializer lists.

In this example, we introduce both the pair-wise :cpp:class:`Lennard-Jones (LJ) model <feasst::ModelLJ>`, and also :cpp:class:`long-range corrections <feasst::LongRangeCorrections>`, which approximately account for the cut off of the LJ potential by assuming a pair-wise radial distance distribution function of unity.
Also note that when creating pointers to FEASST dervied class objects, a convention is to use a helper function which appends the word `Make` onto the class name.

Initialization of the acceptance :cpp:class:`Criteria <feasst::Criteria>` includes thermodynamic constraints such as temperature and chemical potential.

.. literalinclude:: lj_brief.py
   :language: py
   :lines: 9-10
   :lineno-match:

A :cpp:class:`TrialTranslate <feasst::TrialTranslate>` is then introduced which attempts to translate a random particle by a random distance which is bound in each dimension by a `tunable_param`.
This parameter may be adjusted to obtain a desired acceptance ratio, `tunable_target_acceptance`, with the help of :cpp:class:`Tuner <feasst::Tuner>`.

.. literalinclude:: lj_brief.py
   :language: py
   :lines: 11-13
   :lineno-match:

With the help of :cpp:class:`TrialTranslate <feasst::TrialTranslate>`, we can now initialize the number of particles.

.. literalinclude:: lj_brief.py
   :language: py
   :lines: 14
   :lineno-match:

.. warning::

   While :cpp:func:`seek_num_particles <feasst::MonteCarlo::seek_num_particles>` appears simple, care must be taken for when to perform this step.
   A grand canonical simulation is performed here, which is why chemical potential was input for :cpp:class:`Criteria <feasst::Criteria>`.
   This function is likely to change in upcoming versions to reflect more nuanced features which are out of the scope of this tutorial.

Additional :cpp:class:`Analyze <feasst::Analyze>` or :cpp:class:`Modify <feasst::Modify>` may be added at any time to perform some task contingent upon the number of attempted trials.

.. literalinclude:: lj_brief.py
   :language: py
   :lines: 15-20
   :lineno-match:

In this example, :cpp:class:`Log <feasst::Log>` outputs the current status of the trials, :cpp:class:`Movie <feasst::Movie>` outputs the configuration, and :cpp:class:`CheckEnergy <feasst::CheckEnergy>` asserts that the optimized energy calculations match the unoptimized ones.

Finally, a random number generator using the C++ implementation of the Mersenne Twister is seeded by the date.

.. literalinclude:: lj_brief.py
   :language: py
   :lines: 21
   :lineno-match:

.. warning::

   Care must be taken not to run two identical simulations at the same second on an HPC node using the date, or they will have the same seed and may not diverge.
   Instead, consider using a `thread safe` random number generator to seed the simulations.

The simulation is finally run for a number of trial attempts.

.. literalinclude:: lj_brief.py
   :language: py
   :lines: 22
   :lineno-match:

Run the simulation on the command line and you should see something similar to the following::

   [user@host]$ feasst/py/run.sh feasst/tutorial/lj_brief.py
   ...
   energy attempt TrialTranslate tunable
   ...
   -99.7639 980000 0.20339 0.650287
   -127.211 990000 0.186441 0.617773
   -127.339 1000000 0.152542 0.586884

Because the :cpp:class:`Log <feasst::Log>` did not include `file_name` in its :cpp:class:`Arguments <feasst::Arguments>`, it will use standard output to screen.
Thus, you should observe a column of energies, numbers of trials, translation trial acceptance and the value of the tunable parameter.
You should also find the `movie.xyz` trajectory file printed by :cpp:class:`Movie <feasst::Movie>` with an automatically-generated `movie.xyz.vmd` file for use with VMD (e.g., `vmd -e movie.xyz.vmd`).

.. include:: /dev/sphinx/feedback.rst

Additional tutorials
======================

After completing this basic tutorial, check out the tutorials specific to each :doc:`/plugin/README`.
For example, see the tutorials of the :doc:`/plugin/system/README`, :doc:`/plugin/monte_carlo/README` and :doc:`/plugin/flat_histogram/README` plugins.
To complete the tutorials, investigate the contents and comments contained within each script present in the `plugin/<name>/tutorial` directory.
Then run the scripts and see that they produce the expected results.
