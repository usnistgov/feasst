/* This is an interface file for swig.
   This file is created by dev/tools/depend.py . Modifications to this
   file will likely be overwritten in the future. Thus, edit depend.py
   instead.

   usage: /path/to/feasst/py/run.sh /path/to/feasst/dev/tools/depend.py -s /path/to/feasst
 */

%module feasst

%{
#include "configuration/include/properties.h"
#include "configuration/include/typed_entity.h"
#include "configuration/include/bond.h"
#include "monte_carlo/include/tunable.h"
#include "system/include/model.h"
#include "threads/include/thread.h"
#include "threads/include/thread_omp.h"
#include "utils/include/timer.h"
#include "utils/include/utils.h"
#include "utils/include/custom_exception.h"
#include "utils/include/debug.h"
#include "utils/include/utils_file.h"
#include "utils/include/arguments.h"
#include "utils/include/checkpoint.h"
#include "system/include/energy_map.h"
#include "system/include/visit_model_inner.h"
#include "configuration/include/physical_constants.h"
#include "utils/include/progress_report.h"
#include "utils/include/utils_io.h"
#include "utils/include/serialize.h"
#include "utils/include/cache.h"
#include "math/include/table.h"
#include "math/include/accumulator.h"
#include "monte_carlo/include/stepper.h"
#include "math/include/histogram.h"
#include "math/include/formula.h"
#include "shape/include/formula_sine_wave.h"
#include "math/include/formula_polynomial.h"
#include "math/include/utils_math.h"
#include "math/include/position.h"
#include "shape/include/shape.h"
#include "shape/include/slab_sine.h"
#include "shape/include/sphere.h"
#include "shape/include/slab.h"
#include "shape/include/half_space.h"
#include "shape/include/half_space_sine.h"
#include "shape/include/cylinder.h"
#include "shape/include/shape_union.h"
#include "shape/include/shape_intersect.h"
#include "system/include/bond_three_body.h"
#include "system/include/angle_square_well.h"
#include "system/include/bond_two_body.h"
#include "system/include/bond_visitor.h"
#include "system/include/bond_square_well.h"
#include "system/include/neighbor_criteria.h"
#include "cluster/include/energy_map_all.h"
#include "cluster/include/energy_map_neighbor.h"
#include "cluster/include/energy_map_all_criteria.h"
#include "system/include/visit_model.h"
#include "system/include/model_two_body.h"
#include "models/include/yukawa.h"
#include "models/include/square_well.h"
#include "example/include/model_example.h"
#include "system/include/ideal_gas.h"
#include "system/include/model_two_body_factory.h"
#include "system/include/hard_sphere.h"
#include "system/include/lennard_jones.h"
#include "system/include/dont_visit_model.h"
#include "system/include/model_one_body.h"
#include "system/include/model_empty.h"
#include "system/include/model_three_body.h"
#include "system/include/visit_model_intra.h"
#include "system/include/visit_model_cell.h"
#include "system/include/long_range_corrections.h"
#include "system/include/visit_model_bond.h"
#include "configuration/include/site.h"
#include "configuration/include/particle.h"
#include "configuration/include/model_params.h"
#include "system/include/potential.h"
#include "system/include/potential_factory.h"
#include "models/include/lennard_jones_alpha.h"
#include "models/include/lennard_jones_cut_shift.h"
#include "models/include/lennard_jones_force_shift.h"
#include "configuration/include/group.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/select.h"
#include "monte_carlo/include/rosenbluth.h"
#include "monte_carlo/include/acceptance.h"
#include "configuration/include/cells.h"
#include "configuration/include/visit_particles.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/constraint.h"
#include "monte_carlo/include/constrain_num_particles.h"
#include "egce/include/a_equal_b.h"
#include "egce/include/a_half_b.h"
#include "monte_carlo/include/trial_select.h"
#include "monte_carlo/include/trial_select_bond.h"
#include "monte_carlo/include/trial_select_angle.h"
#include "chain/include/select_site_of_type.h"
#include "chain/include/select_perturbed.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "chain/include/select_segment.h"
#include "chain/include/select_end_segment.h"
#include "chain/include/select_reptate.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/select_cluster.h"
#include "cluster/include/select_particle_avb_divalent.h"
#include "monte_carlo/include/perturb.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_remove.h"
#include "monte_carlo/include/trial_factory.h"
#include "cluster/include/trial_avb2.h"
#include "monte_carlo/include/analyze.h"
#include "chain/include/analyze_rigid_bonds.h"
#include "monte_carlo/include/analyze_factory.h"
#include "steppers/include/log_and_movie.h"
#include "steppers/include/wall_clock_limit.h"
#include "steppers/include/check_physicality.h"
#include "steppers/include/check.h"
#include "steppers/include/cpu_time.h"
#include "steppers/include/mean_squared_displacement.h"
#include "steppers/include/num_particles.h"
#include "steppers/include/log.h"
#include "steppers/include/energy.h"
#include "steppers/include/criteria_writer.h"
#include "monte_carlo/include/modify.h"
#include "steppers/include/check_energy.h"
#include "monte_carlo/include/modify_factory.h"
#include "monte_carlo/include/monte_carlo.h"
#include "prefetch/include/prefetch.h"
#include "steppers/include/check_energy_and_tune.h"
#include "steppers/include/tuner.h"
#include "steppers/include/check_properties.h"
#include "steppers/include/criteria_updater.h"
#include "chain/include/recenter_particles.h"
#include "cluster/include/trial_transfer_avb.h"
#include "monte_carlo/include/seek_num_particles.h"
#include "monte_carlo/include/trial_transfer.h"
#include "cluster/include/trial_rigid_cluster.h"
#include "monte_carlo/include/trial_add.h"
#include "cluster/include/trial_avb2_half.h"
#include "cluster/include/trial_remove_avb.h"
#include "cluster/include/trial_add_avb.h"
#include "cluster/include/trial_transfer_avb_divalent.h"
#include "cluster/include/trial_avb4.h"
#include "chain/include/trial_swap_sites.h"
#include "cluster/include/compute_move_cluster.h"
#include "cluster/include/compute_add_avb.h"
#include "cluster/include/compute_remove_avb.h"
#include "cluster/include/compute_gca.h"
#include "cluster/include/compute_add_avb_divalent.h"
#include "cluster/include/compute_remove_avb_divalent.h"
#include "monte_carlo/include/trial_compute_add.h"
#include "monte_carlo/include/trial_compute_remove.h"
#include "monte_carlo/include/trial_compute_move.h"
#include "cluster/include/compute_avb2.h"
#include "cluster/include/compute_avb4.h"
#include "monte_carlo/include/perturb_move.h"
#include "monte_carlo/include/perturb_translate.h"
#include "cluster/include/perturb_point_reflect.h"
#include "monte_carlo/include/perturb_distance.h"
#include "chain/include/perturb_reptate.h"
#include "monte_carlo/include/trial_move.h"
#include "chain/include/trial_pivot.h"
#include "monte_carlo/include/trial_translate.h"
#include "chain/include/trial_reptate.h"
#include "chain/include/trial_crankshaft.h"
#include "system/include/utils.h"
#include "configuration/include/visit_configuration.h"
#include "configuration/include/file_xyz.h"
#include "steppers/include/movie.h"
#include "configuration/include/utils.h"
#include "configuration/include/file_lmp.h"
#include "configuration/include/domain.h"
#include "math/include/random.h"
#include "math/include/constants.h"
#include "math/include/formula_exponential.h"
#include "math/include/matrix.h"
#include "monte_carlo/include/perturb_rotate.h"
#include "chain/include/perturb_pivot.h"
#include "chain/include/perturb_crankshaft.h"
#include "cluster/include/perturb_rotate_com.h"
#include "cluster/include/perturb_move_avb.h"
#include "cluster/include/perturb_add_avb.h"
#include "monte_carlo/include/trial_rotate.h"
#include "monte_carlo/include/perturb_anywhere.h"
#include "chain/include/trial_grow.h"
#include "growth_expanded/include/trial_growth_expanded.h"
#include "chain/include/perturb_site_type.h"
#include "monte_carlo/include/perturb_add.h"
#include "monte_carlo/include/perturb_remove.h"
#include "monte_carlo/include/perturb_distance_angle.h"
#include "math/include/random_mt19937.h"
#include "math/include/solver.h"
#include "math/include/solver_newton_raphson.h"
#include "math/include/solver_bisection.h"
#include "opt_lj/include/visit_model_opt_lj.h"
#include "ewald/include/trial_remove_multiple.h"
#include "ewald/include/trial_add_multiple.h"
#include "ewald/include/charge_self.h"
#include "ewald/include/charge_screened.h"
#include "ewald/include/utils.h"
#include "ewald/include/trial_transfer_multiple.h"
#include "ewald/include/compute_remove_multiple.h"
#include "ewald/include/compute_add_multiple.h"
#include "ewald/include/charge_screened_intra.h"
#include "ewald/include/ewald.h"
#include "ewald/include/check_net_charge.h"
#include "mayer/include/criteria_mayer.h"
#include "mayer/include/trial.h"
#include "confinement/include/model_lj_shape.h"
#include "confinement/include/model_hard_shape.h"
#include "confinement/include/model_table_cartesian.h"
#include "patch/include/patch_angle.h"
#include "patch/include/visit_model_inner_patch.h"
#include "flat_histogram/include/window.h"
#include "flat_histogram/include/window_custom.h"
#include "flat_histogram/include/window_exponential.h"
#include "flat_histogram/include/ln_probability.h"
#include "flat_histogram/include/collection_matrix.h"
#include "flat_histogram/include/bias.h"
#include "flat_histogram/include/transition_matrix.h"
#include "flat_histogram/include/wang_landau.h"
#include "flat_histogram/include/wltm.h"
#include "flat_histogram/include/macrostate.h"
#include "growth_expanded/include/macrostate_growth_expanded.h"
#include "flat_histogram/include/macrostate_num_particles.h"
#include "flat_histogram/include/flat_histogram.h"
#include "flat_histogram/include/clones.h"
using namespace feasst;
%}
%include "std_string.i"
%include "std_vector.i"
%include "std_shared_ptr.i"
%include "std_iostream.i"
%template(IntVector) std::vector<int>;
%template(Int2DVector) std::vector<std::vector<int> >;
%template(DoubleVector) std::vector<double>;
%template(Double2DVector) std::vector<std::vector<double> >;
%template(Double3DVector) std::vector<std::vector<std::vector<double> > >;
using namespace std;
%pythonnondynamic;
%include "std_map.i"
%template(args) std::map<std::string, std::string>;
%template(ModelTwoBodyVector) std::vector<std::shared_ptr<ModelTwoBody> >;
%shared_ptr(feasst::Properties);
%shared_ptr(feasst::PropertiedEntity);
%shared_ptr(feasst::TypedEntity);
%shared_ptr(feasst::Bond);
%shared_ptr(feasst::Angle);
%shared_ptr(feasst::Dihedral);
%shared_ptr(feasst::Improper);
%shared_ptr(feasst::Tunable);
%shared_ptr(feasst::Model);
%shared_ptr(feasst::Thread);
%shared_ptr(feasst::ThreadOMP);
%shared_ptr(feasst::Timer);
%shared_ptr(feasst::CustomException);
%shared_ptr(feasst::Arguments);
%shared_ptr(feasst::Checkpoint);
%shared_ptr(feasst::EnergyMap);
%shared_ptr(feasst::VisitModelInner);
%shared_ptr(feasst::PhysicalConstants);
%shared_ptr(feasst::CODATA2018);
%shared_ptr(feasst::CODATA2014);
%shared_ptr(feasst::CODATA2010);
%shared_ptr(feasst::PhysicalConstantsCustom);
%shared_ptr(feasst::ProgressReport);
%shared_ptr(feasst::Cache);
%shared_ptr(feasst::Table);
%shared_ptr(feasst::Table3D);
%shared_ptr(feasst::Accumulator);
%shared_ptr(feasst::Stepper);
%shared_ptr(feasst::Histogram);
%shared_ptr(feasst::Formula);
%shared_ptr(feasst::FormulaSineWave);
%shared_ptr(feasst::FormulaPolynomial);
%shared_ptr(feasst::Position);
%shared_ptr(feasst::SpatialEntity);
%shared_ptr(feasst::Shape);
%shared_ptr(feasst::ShapedEntity);
%shared_ptr(feasst::SlabSine);
%shared_ptr(feasst::Sphere);
%shared_ptr(feasst::Slab);
%shared_ptr(feasst::HalfSpace);
%shared_ptr(feasst::HalfSpaceSine);
%shared_ptr(feasst::Cylinder);
%shared_ptr(feasst::ShapeUnion);
%shared_ptr(feasst::ShapeIntersect);
%shared_ptr(feasst::BondThreeBody);
%shared_ptr(feasst::AngleSquareWell);
%shared_ptr(feasst::BondTwoBody);
%shared_ptr(feasst::BondVisitor);
%shared_ptr(feasst::BondSquareWell);
%shared_ptr(feasst::NeighborCriteria);
%shared_ptr(feasst::EnergyMapAll);
%shared_ptr(feasst::EnergyMapNeighbor);
%shared_ptr(feasst::EnergyMapAllCriteria);
%shared_ptr(feasst::VisitModel);
%shared_ptr(feasst::ModelTwoBody);
%shared_ptr(feasst::Yukawa);
%shared_ptr(feasst::SquareWell);
%shared_ptr(feasst::ModelExample);
%shared_ptr(feasst::IdealGas);
%shared_ptr(feasst::ModelTwoBodyFactory);
%shared_ptr(feasst::HardSphere);
%shared_ptr(feasst::LennardJones);
%shared_ptr(feasst::DontVisitModel);
%shared_ptr(feasst::ModelOneBody);
%shared_ptr(feasst::ModelEmpty);
%shared_ptr(feasst::ModelThreeBody);
%shared_ptr(feasst::VisitModelIntra);
%shared_ptr(feasst::VisitModelCell);
%shared_ptr(feasst::LongRangeCorrections);
%shared_ptr(feasst::VisitModelBond);
%shared_ptr(feasst::Site);
%shared_ptr(feasst::Particle);
%shared_ptr(feasst::ModelParam);
%shared_ptr(feasst::Epsilon);
%shared_ptr(feasst::Sigma);
%shared_ptr(feasst::CutOff);
%shared_ptr(feasst::Charge);
%shared_ptr(feasst::ModelParams);
%shared_ptr(feasst::Potential);
%shared_ptr(feasst::PotentialFactory);
%shared_ptr(feasst::LennardJonesAlpha);
%shared_ptr(feasst::EnergyAtCutoff);
%shared_ptr(feasst::EnergyDerivAtCutoff);
%shared_ptr(feasst::LennardJonesCutShift);
%shared_ptr(feasst::LennardJonesForceShift);
%shared_ptr(feasst::Group);
%shared_ptr(feasst::ParticleFactory);
%shared_ptr(feasst::Select);
%shared_ptr(feasst::Rosenbluth);
%shared_ptr(feasst::Acceptance);
%shared_ptr(feasst::Cells);
%shared_ptr(feasst::VisitParticles);
%shared_ptr(feasst::LoopOneBody);
%shared_ptr(feasst::Configuration);
%shared_ptr(feasst::System);
%shared_ptr(feasst::Criteria);
%shared_ptr(feasst::Metropolis);
%shared_ptr(feasst::Constraint);
%shared_ptr(feasst::ConstrainNumParticles);
%shared_ptr(feasst::AEqualB);
%shared_ptr(feasst::AHalfB);
%shared_ptr(feasst::TrialSelect);
%shared_ptr(feasst::TrialSelectBond);
%shared_ptr(feasst::TrialSelectAngle);
%shared_ptr(feasst::SelectSiteOfType);
%shared_ptr(feasst::SelectPerturbed);
%shared_ptr(feasst::TrialSelectParticle);
%shared_ptr(feasst::SelectSegment);
%shared_ptr(feasst::SelectEndSegment);
%shared_ptr(feasst::SelectReptate);
%shared_ptr(feasst::SelectParticleAVB);
%shared_ptr(feasst::SelectCluster);
%shared_ptr(feasst::SelectParticleAVBDivalent);
%shared_ptr(feasst::Perturb);
%shared_ptr(feasst::TrialStage);
%shared_ptr(feasst::TrialCompute);
%shared_ptr(feasst::Trial);
%shared_ptr(feasst::TrialRemove);
%shared_ptr(feasst::TrialFactory);
%shared_ptr(feasst::TrialAVB2);
%shared_ptr(feasst::Analyze);
%shared_ptr(feasst::AnalyzeWriteOnly);
%shared_ptr(feasst::AnalyzeUpdateOnly);
%shared_ptr(feasst::AnalyzeRigidBonds);
%shared_ptr(feasst::AnalyzeFactory);
%shared_ptr(feasst::LogAndMovie);
%shared_ptr(feasst::WallClockLimit);
%shared_ptr(feasst::CheckPhysicality);
%shared_ptr(feasst::Check);
%shared_ptr(feasst::CPUTime);
%shared_ptr(feasst::MeanSquaredDisplacement);
%shared_ptr(feasst::NumParticles);
%shared_ptr(feasst::Log);
%shared_ptr(feasst::Energy);
%shared_ptr(feasst::CriteriaWriter);
%shared_ptr(feasst::Modify);
%shared_ptr(feasst::ModifyUpdateOnly);
%shared_ptr(feasst::CheckEnergy);
%shared_ptr(feasst::ModifyFactory);
%shared_ptr(feasst::MonteCarlo);
%shared_ptr(feasst::Pool);
%shared_ptr(feasst::Prefetch);
%shared_ptr(feasst::CheckEnergyAndTune);
%shared_ptr(feasst::Tuner);
%shared_ptr(feasst::CheckProperties);
%shared_ptr(feasst::CriteriaUpdater);
%shared_ptr(feasst::RecenterParticles);
%shared_ptr(feasst::TrialTransferAVB);
%shared_ptr(feasst::SeekNumParticles);
%shared_ptr(feasst::TrialTransfer);
%shared_ptr(feasst::TrialTranslateCluster);
%shared_ptr(feasst::TrialRotateCluster);
%shared_ptr(feasst::TrialRigidCluster);
%shared_ptr(feasst::TrialAdd);
%shared_ptr(feasst::TrialAVB2Half);
%shared_ptr(feasst::TrialRemoveAVB);
%shared_ptr(feasst::TrialAddAVB);
%shared_ptr(feasst::TrialAddAVBDivalent);
%shared_ptr(feasst::TrialRemoveAVBDivalent);
%shared_ptr(feasst::TrialTransferAVBDivalent);
%shared_ptr(feasst::TrialAVB4);
%shared_ptr(feasst::TrialSwapSites);
%shared_ptr(feasst::ComputeMoveCluster);
%shared_ptr(feasst::ComputeAddAVB);
%shared_ptr(feasst::ComputeRemoveAVB);
%shared_ptr(feasst::ComputeGCA);
%shared_ptr(feasst::ComputeAddAVBDivalent);
%shared_ptr(feasst::ComputeRemoveAVBDivalent);
%shared_ptr(feasst::TrialComputeAdd);
%shared_ptr(feasst::TrialComputeRemove);
%shared_ptr(feasst::TrialComputeMove);
%shared_ptr(feasst::ComputeAVB2);
%shared_ptr(feasst::ComputeAVB4);
%shared_ptr(feasst::PerturbMove);
%shared_ptr(feasst::PerturbTranslate);
%shared_ptr(feasst::PerturbPointReflect);
%shared_ptr(feasst::PerturbDistance);
%shared_ptr(feasst::PerturbReptate);
%shared_ptr(feasst::TrialMove);
%shared_ptr(feasst::TrialPivot);
%shared_ptr(feasst::TrialTranslate);
%shared_ptr(feasst::TrialReptate);
%shared_ptr(feasst::TrialCrankshaft);
%shared_ptr(feasst::VisitConfiguration);
%shared_ptr(feasst::LoopConfigOneBody);
%shared_ptr(feasst::FileVMD);
%shared_ptr(feasst::FileXYZ);
%shared_ptr(feasst::Movie);
%shared_ptr(feasst::FileLMP);
%shared_ptr(feasst::Domain);
%shared_ptr(feasst::Random);
%shared_ptr(feasst::FormulaExponential);
%shared_ptr(feasst::Matrix);
%shared_ptr(feasst::MatrixThreeByThree);
%shared_ptr(feasst::RotationMatrix);
%shared_ptr(feasst::PerturbRotate);
%shared_ptr(feasst::PerturbPivot);
%shared_ptr(feasst::PerturbCrankshaft);
%shared_ptr(feasst::PerturbRotateCOM);
%shared_ptr(feasst::PerturbMoveAVB);
%shared_ptr(feasst::PerturbAddAVB);
%shared_ptr(feasst::TrialRotate);
%shared_ptr(feasst::PerturbAnywhere);
%shared_ptr(feasst::TrialGrowLinear);
%shared_ptr(feasst::TrialComputeGrowAdd);
%shared_ptr(feasst::TrialComputeGrowRemove);
%shared_ptr(feasst::TrialComputeGrow);
%shared_ptr(feasst::TrialGrowthExpanded);
%shared_ptr(feasst::PerturbSiteType);
%shared_ptr(feasst::PerturbAdd);
%shared_ptr(feasst::PerturbRemove);
%shared_ptr(feasst::PerturbDistanceAngle);
%shared_ptr(feasst::RandomMT19937);
%shared_ptr(feasst::Solver);
%shared_ptr(feasst::SolverNewtonRaphson);
%shared_ptr(feasst::SolverBisection);
%shared_ptr(feasst::VisitModelOptLJ);
%shared_ptr(feasst::TrialRemoveMultiple);
%shared_ptr(feasst::TrialAddMultiple);
%shared_ptr(feasst::ChargeSelf);
%shared_ptr(feasst::ChargeScreened);
%shared_ptr(feasst::TrialTransferMultiple);
%shared_ptr(feasst::ComputeRemoveMultiple);
%shared_ptr(feasst::ComputeAddMultiple);
%shared_ptr(feasst::ChargeScreenedIntra);
%shared_ptr(feasst::Ewald);
%shared_ptr(feasst::CheckNetCharge);
%shared_ptr(feasst::CriteriaMayer);
%shared_ptr(feasst::TrialComputeMoveMayer);
%shared_ptr(feasst::TrialTranslateMayer);
%shared_ptr(feasst::ModelLJShape);
%shared_ptr(feasst::ModelHardShape);
%shared_ptr(feasst::ModelTableCart3FoldSym);
%shared_ptr(feasst::PatchAngle);
%shared_ptr(feasst::CosPatchAngle);
%shared_ptr(feasst::VisitModelInnerPatch);
%shared_ptr(feasst::Window);
%shared_ptr(feasst::WindowCustom);
%shared_ptr(feasst::WindowExponential);
%shared_ptr(feasst::LnProbability);
%shared_ptr(feasst::TripleBandedCollectionMatrix);
%shared_ptr(feasst::Bias);
%shared_ptr(feasst::TransitionMatrix);
%shared_ptr(feasst::WangLandau);
%shared_ptr(feasst::WLTM);
%shared_ptr(feasst::Macrostate);
%shared_ptr(feasst::MacrostateGrowthExpanded);
%shared_ptr(feasst::MacrostateNumParticles);
%shared_ptr(feasst::FlatHistogram);
%shared_ptr(feasst::Clones);
%include configuration/include/properties.h
%include configuration/include/typed_entity.h
%include configuration/include/bond.h
%include monte_carlo/include/tunable.h
%include system/include/model.h
%include threads/include/thread.h
%include threads/include/thread_omp.h
%include utils/include/timer.h
%include utils/include/utils.h
%include utils/include/custom_exception.h
%include utils/include/debug.h
%include utils/include/utils_file.h
%include utils/include/arguments.h
%include utils/include/checkpoint.h
%include system/include/energy_map.h
%include system/include/visit_model_inner.h
%include configuration/include/physical_constants.h
%include utils/include/progress_report.h
%include utils/include/utils_io.h
%include utils/include/serialize.h
%include utils/include/cache.h
%include math/include/table.h
%include math/include/accumulator.h
%include monte_carlo/include/stepper.h
%include math/include/histogram.h
%include math/include/formula.h
%include shape/include/formula_sine_wave.h
%include math/include/formula_polynomial.h
%include math/include/utils_math.h
%include math/include/position.h
%include shape/include/shape.h
%include shape/include/slab_sine.h
%include shape/include/sphere.h
%include shape/include/slab.h
%include shape/include/half_space.h
%include shape/include/half_space_sine.h
%include shape/include/cylinder.h
%include shape/include/shape_union.h
%include shape/include/shape_intersect.h
%include system/include/bond_three_body.h
%include system/include/angle_square_well.h
%include system/include/bond_two_body.h
%include system/include/bond_visitor.h
%include system/include/bond_square_well.h
%include system/include/neighbor_criteria.h
%include cluster/include/energy_map_all.h
%include cluster/include/energy_map_neighbor.h
%include cluster/include/energy_map_all_criteria.h
%include system/include/visit_model.h
%include system/include/model_two_body.h
%include models/include/yukawa.h
%include models/include/square_well.h
%include example/include/model_example.h
%include system/include/ideal_gas.h
%include system/include/model_two_body_factory.h
%include system/include/hard_sphere.h
%include system/include/lennard_jones.h
%include system/include/dont_visit_model.h
%include system/include/model_one_body.h
%include system/include/model_empty.h
%include system/include/model_three_body.h
%include system/include/visit_model_intra.h
%include system/include/visit_model_cell.h
%include system/include/long_range_corrections.h
%include system/include/visit_model_bond.h
%include configuration/include/site.h
%include configuration/include/particle.h
%include configuration/include/model_params.h
%include system/include/potential.h
%include system/include/potential_factory.h
%include models/include/lennard_jones_alpha.h
%include models/include/lennard_jones_cut_shift.h
%include models/include/lennard_jones_force_shift.h
%include configuration/include/group.h
%include configuration/include/particle_factory.h
%include configuration/include/select.h
%include monte_carlo/include/rosenbluth.h
%include monte_carlo/include/acceptance.h
%include configuration/include/cells.h
%include configuration/include/visit_particles.h
%include configuration/include/configuration.h
%include system/include/system.h
%include monte_carlo/include/criteria.h
%include monte_carlo/include/metropolis.h
%include monte_carlo/include/constraint.h
%include monte_carlo/include/constrain_num_particles.h
%include egce/include/a_equal_b.h
%include egce/include/a_half_b.h
%include monte_carlo/include/trial_select.h
%include monte_carlo/include/trial_select_bond.h
%include monte_carlo/include/trial_select_angle.h
%include chain/include/select_site_of_type.h
%include chain/include/select_perturbed.h
%include monte_carlo/include/trial_select_particle.h
%include chain/include/select_segment.h
%include chain/include/select_end_segment.h
%include chain/include/select_reptate.h
%include cluster/include/select_particle_avb.h
%include cluster/include/select_cluster.h
%include cluster/include/select_particle_avb_divalent.h
%include monte_carlo/include/perturb.h
%include monte_carlo/include/trial_stage.h
%include monte_carlo/include/trial_compute.h
%include monte_carlo/include/trial.h
%include monte_carlo/include/trial_remove.h
%include monte_carlo/include/trial_factory.h
%include cluster/include/trial_avb2.h
%include monte_carlo/include/analyze.h
%include chain/include/analyze_rigid_bonds.h
%include monte_carlo/include/analyze_factory.h
%include steppers/include/log_and_movie.h
%include steppers/include/wall_clock_limit.h
%include steppers/include/check_physicality.h
%include steppers/include/check.h
%include steppers/include/cpu_time.h
%include steppers/include/mean_squared_displacement.h
%include steppers/include/num_particles.h
%include steppers/include/log.h
%include steppers/include/energy.h
%include steppers/include/criteria_writer.h
%include monte_carlo/include/modify.h
%include steppers/include/check_energy.h
%include monte_carlo/include/modify_factory.h
%include monte_carlo/include/monte_carlo.h
%include prefetch/include/prefetch.h
%include steppers/include/check_energy_and_tune.h
%include steppers/include/tuner.h
%include steppers/include/check_properties.h
%include steppers/include/criteria_updater.h
%include chain/include/recenter_particles.h
%include cluster/include/trial_transfer_avb.h
%include monte_carlo/include/seek_num_particles.h
%include monte_carlo/include/trial_transfer.h
%include cluster/include/trial_rigid_cluster.h
%include monte_carlo/include/trial_add.h
%include cluster/include/trial_avb2_half.h
%include cluster/include/trial_remove_avb.h
%include cluster/include/trial_add_avb.h
%include cluster/include/trial_transfer_avb_divalent.h
%include cluster/include/trial_avb4.h
%include chain/include/trial_swap_sites.h
%include cluster/include/compute_move_cluster.h
%include cluster/include/compute_add_avb.h
%include cluster/include/compute_remove_avb.h
%include cluster/include/compute_gca.h
%include cluster/include/compute_add_avb_divalent.h
%include cluster/include/compute_remove_avb_divalent.h
%include monte_carlo/include/trial_compute_add.h
%include monte_carlo/include/trial_compute_remove.h
%include monte_carlo/include/trial_compute_move.h
%include cluster/include/compute_avb2.h
%include cluster/include/compute_avb4.h
%include monte_carlo/include/perturb_move.h
%include monte_carlo/include/perturb_translate.h
%include cluster/include/perturb_point_reflect.h
%include monte_carlo/include/perturb_distance.h
%include chain/include/perturb_reptate.h
%include monte_carlo/include/trial_move.h
%include chain/include/trial_pivot.h
%include monte_carlo/include/trial_translate.h
%include chain/include/trial_reptate.h
%include chain/include/trial_crankshaft.h
%include system/include/utils.h
%include configuration/include/visit_configuration.h
%include configuration/include/file_xyz.h
%include steppers/include/movie.h
%include configuration/include/utils.h
%include configuration/include/file_lmp.h
%include configuration/include/domain.h
%include math/include/random.h
%include math/include/constants.h
%include math/include/formula_exponential.h
%include math/include/matrix.h
%include monte_carlo/include/perturb_rotate.h
%include chain/include/perturb_pivot.h
%include chain/include/perturb_crankshaft.h
%include cluster/include/perturb_rotate_com.h
%include cluster/include/perturb_move_avb.h
%include cluster/include/perturb_add_avb.h
%include monte_carlo/include/trial_rotate.h
%include monte_carlo/include/perturb_anywhere.h
%include chain/include/trial_grow.h
%include growth_expanded/include/trial_growth_expanded.h
%include chain/include/perturb_site_type.h
%include monte_carlo/include/perturb_add.h
%include monte_carlo/include/perturb_remove.h
%include monte_carlo/include/perturb_distance_angle.h
%include math/include/random_mt19937.h
%include math/include/solver.h
%include math/include/solver_newton_raphson.h
%include math/include/solver_bisection.h
%include opt_lj/include/visit_model_opt_lj.h
%include ewald/include/trial_remove_multiple.h
%include ewald/include/trial_add_multiple.h
%include ewald/include/charge_self.h
%include ewald/include/charge_screened.h
%include ewald/include/utils.h
%include ewald/include/trial_transfer_multiple.h
%include ewald/include/compute_remove_multiple.h
%include ewald/include/compute_add_multiple.h
%include ewald/include/charge_screened_intra.h
%include ewald/include/ewald.h
%include ewald/include/check_net_charge.h
%include mayer/include/criteria_mayer.h
%include mayer/include/trial.h
%include confinement/include/model_lj_shape.h
%include confinement/include/model_hard_shape.h
%include confinement/include/model_table_cartesian.h
%include patch/include/patch_angle.h
%include patch/include/visit_model_inner_patch.h
%include flat_histogram/include/window.h
%include flat_histogram/include/window_custom.h
%include flat_histogram/include/window_exponential.h
%include flat_histogram/include/ln_probability.h
%include flat_histogram/include/collection_matrix.h
%include flat_histogram/include/bias.h
%include flat_histogram/include/transition_matrix.h
%include flat_histogram/include/wang_landau.h
%include flat_histogram/include/wltm.h
%include flat_histogram/include/macrostate.h
%include growth_expanded/include/macrostate_growth_expanded.h
%include flat_histogram/include/macrostate_num_particles.h
%include flat_histogram/include/flat_histogram.h
%include flat_histogram/include/clones.h
