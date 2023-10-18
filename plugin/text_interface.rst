***************************
Text Interface
***************************

The first word in each line of a FEASST input text file is the name of a class, followed by pairs of class arguments.
The following related lists of classes and their arguments aid in understanding or modifying :doc:`../tutorial/README`.
The text interface is backwards compatible within major version (e.g., text files using version v0.23.0 should also work with v0.23.1).

.. toctree::

   utils/doc/Checkpoint_arguments

Random Number Generators
=========================

.. toctree::

   math/doc/Random_arguments
   math/doc/RandomMT19937_arguments
   math/doc/RandomModulo_arguments

Configuration
============================

These classes describe the identity of the particles, their positions and the spatial domain in which they reside.

.. toctree::

   configuration/doc/Configuration_arguments
   configuration/doc/Domain_arguments
   configuration/doc/PhysicalConstants
   configuration/doc/NeighborCriteria_arguments

Nonbonded Isotropic Models
======================================

These classes include pair-wise (two-body) isotropic interactions.

.. toctree::

   system/doc/ModelEmpty_arguments
   system/doc/IdealGas_arguments
   system/doc/HardSphere_arguments
   system/doc/LennardJones_arguments
   system/doc/ModelTwoBodyFactory_arguments
   models/doc/SquareWell_arguments
   models/doc/LennardJonesAlpha_arguments
   models/doc/LennardJonesCutShift_arguments
   models/doc/LennardJonesForceShift_arguments
   models/doc/TablePotential_arguments
   models/doc/Mie_arguments
   models/doc/Yukawa_arguments
   charge/doc/Coulomb_arguments
   charge/doc/DebyeHuckel_arguments
   charge/doc/ElectricField_arguments
   confinement/doc/HenryCoefficient_arguments
   example/doc/ModelExample_arguments

Nonbonded Anisotropic Models
=====================================

These classes include anisotropic interactions.

.. toctree::

   patch/doc/VisitModelInnerPatch_arguments
   patch/doc/MoviePatch_arguments
   aniso/doc/VisitModelInnerTable_arguments
   aniso/doc/ContactObjective_arguments

Long-Range Electrostatics
======================================

These classes include Ewald and long-range electrostatics.

.. toctree::

   charge/doc/Ewald_arguments
   charge/doc/ChargeScreened_arguments
   charge/doc/ChargeScreenedIntra_arguments
   charge/doc/ChargeSelf_arguments
   charge/doc/SlabCorrection_arguments

Neighbor lists
===========================

These classes store neighbors and their interaction energies.

.. toctree::

   system/doc/VisitModelInner_arguments
   system/doc/EnergyMap_arguments
   cluster/doc/EnergyMapAll_arguments
   cluster/doc/EnergyMapAllCriteria_arguments
   cluster/doc/EnergyMapNeighbor_arguments
   cluster/doc/EnergyMapNeighborCriteria_arguments

Zero- or One-Body Potentials
==============================

These classes include zero- or one-body interactions.

.. toctree::

   confinement/doc/Background_arguments
   confinement/doc/ModelHardShape_arguments
   confinement/doc/ModelTableCart1DHard_arguments

Nonbonded Potentials
==============================

These classes include various ways that interactions may be computed over a :cpp:class:`Selection <feasst::Select>` or the entire :cpp:class:`Configuration <feasst::Configuration>`.

.. toctree::

   system/doc/Potential_arguments
   system/doc/VisitModel_arguments
   system/doc/DontVisitModel_arguments
   system/doc/VisitModelBond_arguments
   system/doc/VisitModelCell_arguments
   system/doc/VisitModelIntra_arguments
   system/doc/VisitModelIntraMap_arguments
   system/doc/LongRangeCorrections_arguments

Bonded Interactions
=======================

The following bonded interactions are typically specified in fstprt files (see :doc:`../particle/README`):

.. toctree::

   system/doc/RigidBond_arguments
   system/doc/RigidAngle_arguments
   system/doc/RigidDihedral_arguments
   system/doc/BondSquareWell_arguments
   system/doc/AngleSquareWell_arguments
   models/doc/BondHarmonic_arguments
   models/doc/AngleHarmonic_arguments
   models/doc/DihedralHarmonic_arguments
   models/doc/DihedralTraPPE_arguments
   models/doc/FENE_arguments

Thermodynamic Parameters
==========================

.. toctree::

   system/doc/ThermoParams_arguments

Acceptance Criteria
==========================

These classes determine if a :cpp:class:`MonteCarlo <feasst::MonteCarlo>` :cpp:class:`Trial <feasst::Trial>` is accepted or rejected.

.. toctree::

   monte_carlo/doc/Metropolis_arguments
   flat_histogram/doc/FlatHistogram_arguments
   mayer/doc/MayerSampling_arguments
   monte_carlo/doc/AlwaysReject_arguments

Constraints
==========================

These classes impose :cpp:class:`Constraints <feasst::Constraint>` on the acceptance :cpp:class:`Criteria <feasst::Criteria>`.

.. toctree::

   monte_carlo/doc/ConstrainNumParticles_arguments
   egce/doc/AEqualB_arguments
   egce/doc/AHalfB_arguments

Monte Carlo Trials
==========================

These classes add :cpp:class:`Trials <feasst::Trial>`.
:cpp:class:`TrialGrow <feasst::TrialGrow>` has the most features.

.. toctree::

   monte_carlo/doc/TrialTranslate_arguments
   chain/doc/TrialParticlePivot_arguments
   monte_carlo/doc/TrialTransfer_arguments
   chain/doc/TrialGrow_arguments
   monte_carlo/doc/TrialAdd_arguments
   monte_carlo/doc/TrialRemove_arguments
   monte_carlo/doc/TrialVolume_arguments
   beta_expanded/doc/TrialBeta_arguments
   chain/doc/TrialCrankshaft_arguments
   chain/doc/TrialCrankshaftSmall_arguments
   chain/doc/TrialPivot_arguments
   chain/doc/TrialSwapSites_arguments
   charge/doc/TrialAddMultiple_arguments
   charge/doc/TrialRemoveMultiple_arguments
   charge/doc/TrialTransferMultiple_arguments
   cluster/doc/TrialAVB2_arguments
   cluster/doc/TrialAVB4_arguments
   cluster/doc/TrialRigidCluster_arguments
   cluster/doc/TrialTransferAVB_arguments
   cluster/doc/TrialTransferAVBDivalent_arguments
   morph/doc/TrialMorph_arguments

Analyze and Modify
==========================

Analyze happen every so many steps, and do not change the simulation.
Modify happen every so many steps, and may change the simulation.

.. toctree::

   steppers/doc/Check_arguments
   steppers/doc/CheckEnergy_arguments
   steppers/doc/Chirality2D_arguments
   steppers/doc/CPUTime_arguments
   steppers/doc/CriteriaUpdater_arguments
   steppers/doc/CriteriaWriter_arguments
   steppers/doc/Density_arguments
   steppers/doc/DensityProfile_arguments
   steppers/doc/Energy_arguments
   steppers/doc/ExtensiveMoments_arguments
   steppers/doc/Log_arguments
   steppers/doc/MeanSquaredDisplacement_arguments
   steppers/doc/Movie_arguments
   steppers/doc/NumParticles_arguments
   steppers/doc/PairDistributionInner_arguments
   steppers/doc/ProfileTrials_arguments
   steppers/doc/ReadConfigFromFile_arguments
   steppers/doc/Scattering_arguments
   steppers/doc/Tune_arguments
   steppers/doc/Volume_arguments
   steppers/doc/WallClockLimit_arguments
   steppers/doc/WrapParticles_arguments
   chain/doc/AnalyzeBonds_arguments
   chain/doc/EndToEndDistance_arguments
   chain/doc/GhostTrialGrow_arguments
   chain/doc/RadiusOfGyration_arguments
   charge/doc/CheckNetCharge_arguments
   cluster/doc/AnalyzeCluster_arguments
   example/doc/AnalyzeExample_arguments

Actions
==========================

Actions happen one time.

.. toctree::

   monte_carlo/doc/Run_arguments
   monte_carlo/doc/RemoveTrial_arguments
   monte_carlo/doc/RemoveAnalyze_arguments
   monte_carlo/doc/RemoveModify_arguments
   monte_carlo/doc/WriteCheckpoint_arguments
   monte_carlo/doc/WriteModelParams_arguments
   monte_carlo/doc/RefPotential_arguments
   monte_carlo/doc/ConvertToRefPotential_arguments
   example/doc/ActionExample_arguments

Flat Histogram
==========================

Flat-histogram simulations bias along a 1-dimensional macrostate.

.. toctree::

   flat_histogram/doc/FlatHistogram_arguments
   flat_histogram/doc/MacrostateNumParticles_arguments
   flat_histogram/doc/MacrostateEnergy_arguments
   beta_expanded/doc/MacrostateBeta_arguments
   flat_histogram/doc/WangLandau_arguments
   flat_histogram/doc/TransitionMatrix_arguments
   flat_histogram/doc/CollectionMatrix_arguments
   flat_histogram/doc/WLTM_arguments
   flat_histogram/doc/CollectionMatrixSplice_arguments
   flat_histogram/doc/Window_arguments
   flat_histogram/doc/WindowExponential_arguments
   flat_histogram/doc/WindowCustom_arguments

Gibbs Ensemble
===================

Simultaneous simulate multiple configurations and transfer particles and/or volume between them.

.. toctree::

   gibbs/doc/TrialGibbsVolumeTransfer_arguments
   gibbs/doc/TrialGibbsParticleTransfer_arguments
   gibbs/doc/PressureFromTestVolume_arguments
   gibbs/doc/CheckConstantVolume_arguments

The information available here is identical to the :doc:`README` documentation, but with only the class arguments and no member functions.
The list above is pruned manually from the :doc:`README` documentation.
Therefore the :doc:`README` documentation represents all capabilities of FEASST, some of which are hidden from text users for various reasons (e.g., because the classes may be used internally, still in development, buggy, depreciated, debugging tools or simply forgotten).

