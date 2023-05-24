***********
Monte Carlo
***********

Monte Carlo randomly samples the system using trials which may be accepted or rejected according to the acceptance criteria.

.. figure:: doc/monte_carlo.svg

:cpp:class:`MonteCarlo <feasst::MonteCarlo>` contains :cpp:class:`System <feasst::System>`, acceptance :cpp:class:`Criteria <feasst::Criteria>`, :cpp:class:`Trials <feasst::Trial>` and, finally, the usually infrequent :cpp:class:`Analyze <feasst::Analyze>` and :cpp:class:`Modify <feasst::Modify>` between :cpp:class:`Trials <feasst::Trial>`.
:cpp:class:`System <feasst::System>` contains :cpp:class:`Configurations <feasst::Configuration>` and :cpp:class:`Potentials <feasst::Potential>`.
:cpp:class:`Configurations <feasst::Configuration>` contain :cpp:class:`Particles <feasst::Particle>` and the :cpp:class:`Domain <feasst::Doman>`.
The interaction :cpp:class:`Potentials <feasst::Potential>` contain :cpp:class:`Models <feasst::Model>` and :cpp:class:`VisitModels <feasst::VisitModel>`.
The acceptance :cpp:class:`Criteria <feasst::Criteria>` may be either :cpp:class:`Metropolis <feasst::Metropolis>` or :cpp:class:`FlatHistogram <feasst::FlatHistogram>`.

A :cpp:class:`Trial <feasst::Trial>` contains the following objects.

.. figure:: doc/trial.svg

Why even have molecular simulations in the first place? Inspired by the second figure in Allen and Tildesley:

.. figure:: doc/sim.svg

Tutorial
=========

.. toctree::
   :glob:

   tutorial/tutorial*

FEASST plugin dependencies
============================

* :doc:`../system/README`

API
===

.. toctree::

   doc/toc

