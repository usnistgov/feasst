***************************
Text Interface
***************************

The syntax of a FEASST text input file is follows.
Each line begins with the name of a class in FEASST, listed below.
The class name is then proceeded by space-separated pairs of argument names and values.
Minor version changes support the same text files, but major version changes do not.
For example, all text files for version v0.23.0 should also work with v0.23.1, but may not work with v0.24.0.
For updating to a newer version, see the change log below.

.. toctree::
   :maxdepth: 1

   utils/doc/Checkpoint_arguments

Random Number Generators
=========================

These are the number generators available for use.

.. toctree::
   :maxdepth: 1

   math/doc/RandomMT19937_arguments
   math/doc/RandomModulo_arguments

Configuration
============================

These classes describe the identity of the particles, their positions and the spatial domain in which they reside.

.. toctree::
   :maxdepth: 1

   configuration/doc/Configuration_arguments
   configuration/doc/Domain_arguments
   configuration/doc/PhysicalConstants

Nonbonded Isotropic Models
======================================

These classes include pair-wise (two-body) isotropic interactions models.

.. toctree::
   :maxdepth: 1

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
   models/doc/Jagla_arguments
   charge/doc/Coulomb_arguments
   charge/doc/DebyeHuckel_arguments
   charge/doc/ElectricField_arguments
   confinement/doc/HenryCoefficient_arguments
   server/doc/ModelServer_arguments
   example/doc/ModelExample_arguments

Nonbonded Potentials
==============================

These classes include various ways that interactions may be computed over a :cpp:class:`Selection <feasst::Select>` or the entire :cpp:class:`Configuration <feasst::Configuration>`.

.. toctree::
   :maxdepth: 1

   system/doc/Potential_arguments
   system/doc/VisitModel_arguments
   system/doc/DontVisitModel_arguments
   system/doc/VisitModelBond_arguments
   system/doc/VisitModelCell_arguments
   system/doc/VisitModelIntra_arguments
   system/doc/VisitModelIntraMap_arguments
   system/doc/LongRangeCorrections_arguments

Nonbonded Anisotropic Models
=====================================

These classes include anisotropic interactions.

.. toctree::
   :maxdepth: 1

   patch/doc/VisitModelInnerPatch_arguments
   patch/doc/MoviePatch_arguments
   aniso/doc/VisitModelInnerTable_arguments
   aniso/doc/TabulateTwoRigidBody3D_arguments

Long-Range Electrostatics
======================================

These classes include Ewald and long-range electrostatics.

.. toctree::
   :maxdepth: 1

   charge/doc/Ewald_arguments
   charge/doc/ChargeScreened_arguments
   charge/doc/ChargeScreenedIntra_arguments
   charge/doc/ChargeSelf_arguments
   charge/doc/SlabCorrection_arguments

Neighbor Lists
===========================

These classes store neighbors and their interaction energies.

.. toctree::
   :maxdepth: 1

   configuration/doc/NeighborCriteria_arguments
   system/doc/VisitModelInner_arguments
   system/doc/EnergyMap_arguments
   cluster/doc/EnergyMapAll_arguments
   cluster/doc/EnergyMapAllCriteria_arguments
   cluster/doc/EnergyMapNeighbor_arguments
   cluster/doc/EnergyMapNeighborCriteria_arguments

One-Body Potentials
==============================

These classes include zero- and one-body interactions.

.. toctree::
   :maxdepth: 1

   confinement/doc/Background_arguments
   confinement/doc/ModelHardShape_arguments
   confinement/doc/ModelLJShape_arguments
   confinement/doc/ModelTableCart1DHard_arguments
   shape/doc/ShapeFile_arguments
   shape/doc/Cuboid_arguments
   shape/doc/Cylinder_arguments
   shape/doc/FiniteCylinder_arguments
   shape/doc/FormulaSineWave_arguments
   shape/doc/HalfSpace_arguments
   shape/doc/Slab_arguments
   shape/doc/SlabSine_arguments
   shape/doc/Sphere_arguments

Bonded Interactions
=======================

The following bonded interactions are typically specified in fstprt files (see :doc:`../particle/README`):

.. toctree::
   :maxdepth: 1

   system/doc/RigidBond_arguments
   system/doc/RigidAngle_arguments
   system/doc/RigidDihedral_arguments
   system/doc/BondSquareWell_arguments
   system/doc/AngleSquareWell_arguments
   models/doc/BondHarmonic_arguments
   models/doc/AngleHarmonic_arguments
   models/doc/DihedralHarmonic_arguments
   models/doc/DihedralTraPPE_arguments
   models/doc/DihedralRyckaertBellemans_arguments
   models/doc/FENE_arguments

Thermodynamic Parameters
==========================

.. toctree::
   :maxdepth: 1

   system/doc/ThermoParams_arguments

Acceptance Criteria
==========================

These classes determine if a :cpp:class:`MonteCarlo <feasst::MonteCarlo>` :cpp:class:`Trial <feasst::Trial>` is accepted or rejected.

.. toctree::
   :maxdepth: 1

   monte_carlo/doc/Metropolis_arguments
   flat_histogram/doc/FlatHistogram_arguments
   mayer/doc/MayerSampling_arguments
   monte_carlo/doc/AlwaysReject_arguments

Constraints
==========================

These classes impose :cpp:class:`Constraints <feasst::Constraint>` on the acceptance :cpp:class:`Criteria <feasst::Criteria>`.

.. toctree::
   :maxdepth: 1

   monte_carlo/doc/ConstrainNumParticles_arguments
   egce/doc/AEqualB_arguments
   egce/doc/AHalfB_arguments

Monte Carlo Trials
==========================

These classes add :cpp:class:`Trials <feasst::Trial>`.

.. toctree::
   :maxdepth: 1

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
   cluster/doc/TrialTranslateCluster_arguments
   cluster/doc/TrialRotateCluster_arguments
   cluster/doc/TrialTransferAVB_arguments
   cluster/doc/TrialAddAVB_arguments
   cluster/doc/TrialRemoveAVB_arguments
   cluster/doc/TrialTransferAVBDivalent_arguments
   cluster/doc/TrialAddAVBDivalent_arguments
   cluster/doc/TrialRemoveAVBDivalent_arguments
   morph/doc/TrialMorph_arguments

Analyze
==========================

:cpp:class:`Analyze <feasst::Analyze>` update/write every fixed number of trials, and do not change the simulation.

.. toctree::
   :maxdepth: 1

   steppers/doc/Check_arguments
   steppers/doc/Chirality2D_arguments
   steppers/doc/CPUTime_arguments
   steppers/doc/CriteriaWriter_arguments
   steppers/doc/Density_arguments
   steppers/doc/DensityProfile_arguments
   steppers/doc/Energy_arguments
   steppers/doc/ExtensiveMoments_arguments
   steppers/doc/HeatCapacity_arguments
   steppers/doc/Log_arguments
   steppers/doc/MeanSquaredDisplacement_arguments
   steppers/doc/Movie_arguments
   steppers/doc/NumParticles_arguments
   steppers/doc/ProfileTrials_arguments
   steppers/doc/Scattering_arguments
   steppers/doc/Volume_arguments
   steppers/doc/WallClockLimit_arguments
   chain/doc/AnalyzeBonds_arguments
   chain/doc/EndToEndDistance_arguments
   chain/doc/RadiusOfGyration_arguments
   charge/doc/CheckNetCharge_arguments
   cluster/doc/AnalyzeCluster_arguments
   example/doc/AnalyzeExample_arguments

Modify
============================

:cpp:class:`Modify <feasst::Modify>` update/write every fixed number of trials, but may change the simulation.

.. toctree::
   :maxdepth: 1

   steppers/doc/CheckEnergy_arguments
   steppers/doc/CriteriaUpdater_arguments
   steppers/doc/PairDistributionInner_arguments
   steppers/doc/ReadConfigFromFile_arguments
   steppers/doc/Tune_arguments
   steppers/doc/WrapParticles_arguments
   chain/doc/GhostTrialGrow_arguments

Actions
==========================

:cpp:class:`Action <feasst::Action>` s happen once.

.. toctree::
   :maxdepth: 1

   monte_carlo/doc/Run_arguments
   monte_carlo/doc/RemoveTrial_arguments
   monte_carlo/doc/RemoveAnalyze_arguments
   monte_carlo/doc/RemoveModify_arguments
   monte_carlo/doc/WriteCheckpoint_arguments
   monte_carlo/doc/WriteModelParams_arguments
   monte_carlo/doc/RefPotential_arguments
   monte_carlo/doc/ConvertToRefPotential_arguments
   steppers/doc/WriteStepper_arguments
   example/doc/ActionExample_arguments

Flat Histogram
==========================

:cpp:class:`FlatHistogram <feasst::FlatHistogram>` simulations :cpp:class:`Bias <feasst::Bias>` along a 1-dimensional :cpp:class:`Macrostate <feasst::Macrostate>`.

.. toctree::
   :maxdepth: 1

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

Simulate multiple :cpp:class:`Configuration <feasst::Configuration>` s and transfer particles and/or volume between them.

.. toctree::
   :maxdepth: 1

   gibbs/doc/TrialGibbsVolumeTransfer_arguments
   gibbs/doc/TrialGibbsParticleTransfer_arguments
   gibbs/doc/PressureFromTestVolume_arguments
   gibbs/doc/CheckConstantVolume_arguments

Change Log
================

Below is a list of all notable changes made to the text interface which will likely lead to the errors if older text interface scripts use the newer version.
Renamed arguments are shown as Class::old_argument->new_argument.

v0.24.6 to v0.25.0
----------------------

* Mie::n,m arguments replaced with mie_lambda_r,mie_lambda_a "Site Properties"
* TablePotential base class VisitModel->Model
* Trial::weight_per_number->weight_per_number_fraction

v0.23.1 to v0.24.0
----------------------

* Stepper::file_name->output_file
* WriteModelParams::file_name->output_file
* Checkpoint::file_name->checkpoint_file
* GhostTrialGrow::trial_grow_file->grow_file
* TrialGrowFile::file_name->grow_file
* ReadConfigFromFile::file_name->input_file
* ShapeFile::file_name->shape_file
* ModelTableCart3DIntegr::file_name->table_file
* ModelTableCart3DIntegr::shape_file_name->shape_file

The information available above is pruned from the :doc:`README` documentation, with only the class arguments and no member functions.
Therefore the :doc:`README` documentation represents all capabilities of FEASST, many of which are hidden from text users for various reasons (e.g., because the classes may be used only internally, in development, buggy, deprecated or debugging tools).

