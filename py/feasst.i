/* This is an interface file for swig.
   This file is created by dev/tools/depend.py . Modifications to this
   file will likely be overwritten in the future. Thus, edit depend.py
   instead.

   usage: /path/to/feasst/py/run.sh /path/to/feasst/dev/tools/depend.py -s /path/to/feasst
 */

%module(directors="1") feasst

%feature("director:except") {
  if( $error != NULL ) {
    PyObject *ptype, *pvalue, *ptraceback;
    PyErr_Fetch( &ptype, &pvalue, &ptraceback );
    PyErr_Restore( ptype, pvalue, ptraceback );
    PyErr_Print();
    Py_Exit(1);
  }
}

%{
#include "configuration/include/properties.h"
#include "configuration/include/typed_entity.h"
#include "configuration/include/bond.h"
#include "monte_carlo/include/tunable.h"
#include "system/include/model.h"
#include "system/include/synchronize_data.h"
#include "utils/include/timer.h"
#include "utils/include/file.h"
#include "utils/include/argument_parse.h"
#include "utils/include/io.h"
#include "utils/include/utils.h"
#include "utils/include/custom_exception.h"
#include "utils/include/debug.h"
#include "utils/include/arguments.h"
#include "utils/include/checkpoint.h"
#include "utils/include/serialize.h"
#include "system/include/energy_map.h"
#include "system/include/visit_model_inner.h"
#include "system/include/thermo_params.h"
#include "configuration/include/physical_constants.h"
#include "utils/include/progress_report.h"
#include "utils/include/cache.h"
#include "math/include/table.h"
#include "math/include/accumulator.h"
#include "monte_carlo/include/stepper.h"
#include "math/include/histogram.h"
#include "math/include/formula.h"
#include "math/include/formula_polynomial.h"
#include "math/include/utils_math.h"
#include "math/include/position.h"
#include "system/include/bond_three_body.h"
#include "system/include/angle_square_well.h"
#include "system/include/bond_two_body.h"
#include "system/include/bond_visitor.h"
#include "system/include/bond_square_well.h"
#include "system/include/neighbor_criteria.h"
#include "cluster/include/energy_map_all.h"
#include "cluster/include/energy_map_neighbor.h"
#include "cluster/include/energy_map_neighbor_criteria.h"
#include "cluster/include/energy_map_all_criteria.h"
#include "system/include/visit_model.h"
#include "system/include/model_two_body.h"
#include "models/include/yukawa.h"
#include "models/include/square_well.h"
#include "system/include/ideal_gas.h"
#include "system/include/model_two_body_factory.h"
#include "system/include/hard_sphere.h"
#include "system/include/lennard_jones.h"
#include "system/include/dont_visit_model.h"
#include "system/include/model_one_body.h"
#include "system/include/model_empty.h"
#include "system/include/model_three_body.h"
#include "system/include/visit_model_intra_map.h"
#include "system/include/visit_model_intra.h"
#include "system/include/visit_model_cell.h"
#include "system/include/long_range_corrections.h"
#include "system/include/visit_model_bond.h"
#include "configuration/include/site.h"
#include "configuration/include/particle.h"
#include "configuration/include/model_params.h"
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
#include "system/include/potential.h"
#include "system/include/potential_factory.h"
#include "system/include/system.h"
#include "system/include/utils.h"
#include "monte_carlo/include/trial_select.h"
#include "monte_carlo/include/trial_select_bond.h"
#include "monte_carlo/include/trial_select_angle.h"
#include "monte_carlo/include/perturb.h"
#include "monte_carlo/include/perturb_move.h"
#include "monte_carlo/include/perturb_translate.h"
#include "monte_carlo/include/perturb_distance.h"
#include "cluster/include/perturb_point_reflect.h"
#include "monte_carlo/include/perturb_volume.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/select_particle_avb_divalent.h"
#include "cluster/include/select_cluster.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute.h"
#include "monte_carlo/include/trial.h"
#include "cluster/include/trial_avb4.h"
#include "monte_carlo/include/trial_factory.h"
#include "cluster/include/trial_avb2.h"
#include "cluster/include/trial_transfer_avb_divalent.h"
#include "cluster/include/trial_transfer_avb.h"
#include "monte_carlo/include/analyze.h"
#include "monte_carlo/include/analyze_factory.h"
#include "steppers/include/log_and_movie.h"
#include "steppers/include/seek_analyze.h"
#include "steppers/include/wall_clock_limit.h"
#include "steppers/include/profile_trials.h"
#include "steppers/include/volume.h"
#include "steppers/include/check_physicality.h"
#include "steppers/include/check.h"
#include "steppers/include/extensive_moments.h"
#include "steppers/include/cpu_time.h"
#include "steppers/include/mean_squared_displacement.h"
#include "steppers/include/num_particles.h"
#include "steppers/include/log.h"
#include "steppers/include/energy.h"
#include "steppers/include/criteria_writer.h"
#include "monte_carlo/include/modify.h"
#include "steppers/include/pair_distribution.h"
#include "steppers/include/check_energy.h"
#include "monte_carlo/include/modify_factory.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/check_energy_and_tune.h"
#include "steppers/include/seek_modify.h"
#include "steppers/include/tuner.h"
#include "steppers/include/check_properties.h"
#include "steppers/include/criteria_updater.h"
#include "monte_carlo/include/seek_num_particles.h"
#include "monte_carlo/include/trials.h"
#include "cluster/include/trial_rigid_cluster.h"
#include "monte_carlo/include/trial_compute_volume.h"
#include "monte_carlo/include/trial_compute_remove.h"
#include "cluster/include/compute_move_cluster.h"
#include "cluster/include/compute_add_avb.h"
#include "cluster/include/compute_remove_avb.h"
#include "cluster/include/compute_gca.h"
#include "monte_carlo/include/trial_compute_add.h"
#include "cluster/include/compute_add_avb_divalent.h"
#include "monte_carlo/include/trial_compute_move.h"
#include "monte_carlo/include/trial_move.h"
#include "cluster/include/compute_remove_avb_divalent.h"
#include "cluster/include/compute_avb2.h"
#include "cluster/include/compute_avb4.h"
#include "monte_carlo/include/constraint.h"
#include "monte_carlo/include/constrain_num_particles.h"
#include "configuration/include/visit_configuration.h"
#include "configuration/include/file_xyz.h"
#include "steppers/include/movie.h"
#include "configuration/include/utils.h"
#include "configuration/include/file_lmp.h"
#include "configuration/include/domain.h"
#include "math/include/random.h"
#include "math/include/constants.h"
#include "math/include/minimize.h"
#include "math/include/golden_search.h"
#include "math/include/quadratic_equation.h"
#include "math/include/formula_exponential.h"
#include "math/include/matrix.h"
#include "monte_carlo/include/perturb_rotate.h"
#include "cluster/include/perturb_rotate_com.h"
#include "cluster/include/perturb_move_avb.h"
#include "cluster/include/perturb_add_avb.h"
#include "monte_carlo/include/perturb_anywhere.h"
#include "monte_carlo/include/perturb_add.h"
#include "monte_carlo/include/perturb_remove.h"
#include "monte_carlo/include/perturb_distance_angle.h"
#include "math/include/random_mt19937.h"
#include "math/include/solver.h"
#include "math/include/solver_newton_raphson.h"
#include "math/include/solver_bisection.h"
#include "math/include/solver_brent_dekker.h"
using namespace feasst;
%}
%include "std_string.i"
%include "std_vector.i"
%include "std_shared_ptr.i"
%include "std_iostream.i"
%include "stdint.i"
%template(IntVector) std::vector<int>;
%template(Int2DVector) std::vector<std::vector<int> >;
%template(DoubleVector) std::vector<double>;
%template(Double2DVector) std::vector<std::vector<double> >;
%template(Double3DVector) std::vector<std::vector<std::vector<double> > >;
using namespace std;
%pythonnondynamic;
%include "std_map.i"
%template(args) std::map<std::string, std::string>;
%template(ArgsVector) std::vector<std::map<std::string, std::string> >;
%template(ModelTwoBodyVector) std::vector<std::shared_ptr<ModelTwoBody> >;
%feature("director") Potential;
%include "std_pair.i"
%template(Map2) std::vector<std::pair<int, std::vector<double>>>;
%template(Map3) std::vector<std::pair<int, std::vector<std::pair<int, std::vector<double>>>>>;
%template(Map4) std::vector<std::pair<int, std::vector<std::pair<int, std::vector<std::pair<int, std::vector<double>>>>>>>;
%template(MapNew) std::vector<std::pair<int, std::vector<std::pair<int, std::vector<std::pair<int, std::vector<std::pair<int, std::vector<double>>>>>>>>>;
%template(MapOld) std::vector<std::vector<std::vector<std::pair<int, std::vector<std::pair<int, std::vector<double>>>>>>>;
%shared_ptr(feasst::Properties);
%shared_ptr(feasst::PropertiedEntity);
%shared_ptr(feasst::TypedEntity);
%shared_ptr(feasst::Bond);
%shared_ptr(feasst::Angle);
%shared_ptr(feasst::Dihedral);
%shared_ptr(feasst::Improper);
%shared_ptr(feasst::Tunable);
%shared_ptr(feasst::Model);
%shared_ptr(feasst::SynchronizeData);
%shared_ptr(feasst::Timer);
%shared_ptr(feasst::ArgumentParse);
%shared_ptr(feasst::CustomException);
%shared_ptr(feasst::Arguments);
%shared_ptr(feasst::Checkpoint);
%shared_ptr(feasst::EnergyMap);
%shared_ptr(feasst::VisitModelInner);
%shared_ptr(feasst::ThermoParams);
%shared_ptr(feasst::PhysicalConstants);
%shared_ptr(feasst::CODATA2018);
%shared_ptr(feasst::CODATA2014);
%shared_ptr(feasst::CODATA2010);
%shared_ptr(feasst::PhysicalConstantsCustom);
%shared_ptr(feasst::ProgressReport);
%shared_ptr(feasst::Cache);
%shared_ptr(feasst::Table);
%shared_ptr(feasst::Table1D);
%shared_ptr(feasst::Table2D);
%shared_ptr(feasst::Table3D);
%shared_ptr(feasst::Accumulator);
%shared_ptr(feasst::Stepper);
%shared_ptr(feasst::Histogram);
%shared_ptr(feasst::Formula);
%shared_ptr(feasst::FormulaPolynomial);
%shared_ptr(feasst::Position);
%shared_ptr(feasst::SpatialEntity);
%shared_ptr(feasst::BondThreeBody);
%shared_ptr(feasst::AngleSquareWell);
%shared_ptr(feasst::BondTwoBody);
%shared_ptr(feasst::BondVisitor);
%shared_ptr(feasst::BondSquareWell);
%shared_ptr(feasst::NeighborCriteria);
%shared_ptr(feasst::EnergyMapAll);
%shared_ptr(feasst::EnergyMapNeighbor);
%shared_ptr(feasst::EnergyMapNeighborCriteria);
%shared_ptr(feasst::EnergyMapAllCriteria);
%shared_ptr(feasst::VisitModel);
%shared_ptr(feasst::ModelTwoBody);
%shared_ptr(feasst::Yukawa);
%shared_ptr(feasst::SquareWell);
%shared_ptr(feasst::IdealGas);
%shared_ptr(feasst::ModelTwoBodyFactory);
%shared_ptr(feasst::HardSphere);
%shared_ptr(feasst::LennardJones);
%shared_ptr(feasst::DontVisitModel);
%shared_ptr(feasst::ModelOneBody);
%shared_ptr(feasst::ModelEmpty);
%shared_ptr(feasst::ModelThreeBody);
%shared_ptr(feasst::VisitModelIntraMap);
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
%shared_ptr(feasst::Potential);
%shared_ptr(feasst::PotentialFactory);
%shared_ptr(feasst::System);
%shared_ptr(feasst::TrialSelect);
%shared_ptr(feasst::TrialSelectBond);
%shared_ptr(feasst::TrialSelectAngle);
%shared_ptr(feasst::Perturb);
%shared_ptr(feasst::PerturbMove);
%shared_ptr(feasst::PerturbTranslate);
%shared_ptr(feasst::PerturbDistance);
%shared_ptr(feasst::PerturbPointReflect);
%shared_ptr(feasst::PerturbVolume);
%shared_ptr(feasst::TrialSelectParticle);
%shared_ptr(feasst::SelectParticleAVB);
%shared_ptr(feasst::SelectParticleAVBDivalent);
%shared_ptr(feasst::SelectCluster);
%shared_ptr(feasst::Criteria);
%shared_ptr(feasst::Metropolis);
%shared_ptr(feasst::TrialStage);
%shared_ptr(feasst::TrialCompute);
%shared_ptr(feasst::Trial);
%shared_ptr(feasst::TrialFactory);
%shared_ptr(feasst::Analyze);
%shared_ptr(feasst::AnalyzeWriteOnly);
%shared_ptr(feasst::AnalyzeUpdateOnly);
%shared_ptr(feasst::AnalyzeFactory);
%shared_ptr(feasst::LogAndMovie);
%shared_ptr(feasst::AnalyzeData);
%shared_ptr(feasst::AccumulatorAverage);
%shared_ptr(feasst::AccumulatorSum);
%shared_ptr(feasst::AccumulatorSumOfSquared);
%shared_ptr(feasst::AccumulatorMoment);
%shared_ptr(feasst::SeekAnalyze);
%shared_ptr(feasst::WallClockLimit);
%shared_ptr(feasst::ProfileTrials);
%shared_ptr(feasst::Volume);
%shared_ptr(feasst::CheckPhysicality);
%shared_ptr(feasst::Check);
%shared_ptr(feasst::ExtensiveMoments);
%shared_ptr(feasst::CPUTime);
%shared_ptr(feasst::MeanSquaredDisplacement);
%shared_ptr(feasst::NumParticles);
%shared_ptr(feasst::Log);
%shared_ptr(feasst::Energy);
%shared_ptr(feasst::CriteriaWriter);
%shared_ptr(feasst::Modify);
%shared_ptr(feasst::ModifyUpdateOnly);
%shared_ptr(feasst::PairDistributionInner);
%shared_ptr(feasst::PairDistribution);
%shared_ptr(feasst::CheckEnergy);
%shared_ptr(feasst::ModifyFactory);
%shared_ptr(feasst::MonteCarlo);
%shared_ptr(feasst::CheckEnergyAndTune);
%shared_ptr(feasst::SeekModify);
%shared_ptr(feasst::Tuner);
%shared_ptr(feasst::CheckProperties);
%shared_ptr(feasst::CriteriaUpdater);
%shared_ptr(feasst::SeekNumParticles);
%shared_ptr(feasst::TrialComputeVolume);
%shared_ptr(feasst::TrialComputeRemove);
%shared_ptr(feasst::ComputeMoveCluster);
%shared_ptr(feasst::ComputeAddAVB);
%shared_ptr(feasst::ComputeRemoveAVB);
%shared_ptr(feasst::ComputeGCA);
%shared_ptr(feasst::TrialComputeAdd);
%shared_ptr(feasst::ComputeAddAVBDivalent);
%shared_ptr(feasst::TrialComputeMove);
%shared_ptr(feasst::ComputeRemoveAVBDivalent);
%shared_ptr(feasst::ComputeAVB2);
%shared_ptr(feasst::ComputeAVB4);
%shared_ptr(feasst::Constraint);
%shared_ptr(feasst::ConstrainNumParticles);
%shared_ptr(feasst::VisitConfiguration);
%shared_ptr(feasst::LoopConfigOneBody);
%shared_ptr(feasst::FileVMD);
%shared_ptr(feasst::FileXYZ);
%shared_ptr(feasst::Movie);
%shared_ptr(feasst::FileLMP);
%shared_ptr(feasst::Domain);
%shared_ptr(feasst::Random);
%shared_ptr(feasst::Minimize);
%shared_ptr(feasst::GoldenSearch);
%shared_ptr(feasst::FormulaExponential);
%shared_ptr(feasst::Matrix);
%shared_ptr(feasst::RotationMatrix);
%shared_ptr(feasst::PerturbRotate);
%shared_ptr(feasst::PerturbRotateCOM);
%shared_ptr(feasst::PerturbMoveAVB);
%shared_ptr(feasst::PerturbAddAVB);
%shared_ptr(feasst::PerturbAnywhere);
%shared_ptr(feasst::PerturbAdd);
%shared_ptr(feasst::PerturbRemove);
%shared_ptr(feasst::PerturbDistanceAngle);
%shared_ptr(feasst::RandomMT19937);
%shared_ptr(feasst::Solver);
%shared_ptr(feasst::SolverNewtonRaphson);
%shared_ptr(feasst::SolverBisection);
%shared_ptr(feasst::SolverBrentDekker);
%include configuration/include/properties.h
%include configuration/include/typed_entity.h
%include configuration/include/bond.h
%include monte_carlo/include/tunable.h
%include system/include/model.h
%include system/include/synchronize_data.h
%include utils/include/timer.h
%include utils/include/file.h
%include utils/include/argument_parse.h
%include utils/include/io.h
%include utils/include/utils.h
%include utils/include/custom_exception.h
%include utils/include/debug.h
%include utils/include/arguments.h
%include utils/include/checkpoint.h
%include utils/include/serialize.h
%include system/include/energy_map.h
%include system/include/visit_model_inner.h
%include system/include/thermo_params.h
%include configuration/include/physical_constants.h
%include utils/include/progress_report.h
%include utils/include/cache.h
%include math/include/table.h
%include math/include/accumulator.h
%include monte_carlo/include/stepper.h
%include math/include/histogram.h
%include math/include/formula.h
%include math/include/formula_polynomial.h
%include math/include/utils_math.h
%include math/include/position.h
%include system/include/bond_three_body.h
%include system/include/angle_square_well.h
%include system/include/bond_two_body.h
%include system/include/bond_visitor.h
%include system/include/bond_square_well.h
%include system/include/neighbor_criteria.h
%include cluster/include/energy_map_all.h
%include cluster/include/energy_map_neighbor.h
%include cluster/include/energy_map_neighbor_criteria.h
%include cluster/include/energy_map_all_criteria.h
%include system/include/visit_model.h
%include system/include/model_two_body.h
%include models/include/yukawa.h
%include models/include/square_well.h
%include system/include/ideal_gas.h
%include system/include/model_two_body_factory.h
%include system/include/hard_sphere.h
%include system/include/lennard_jones.h
%include system/include/dont_visit_model.h
%include system/include/model_one_body.h
%include system/include/model_empty.h
%include system/include/model_three_body.h
%include system/include/visit_model_intra_map.h
%include system/include/visit_model_intra.h
%include system/include/visit_model_cell.h
%include system/include/long_range_corrections.h
%include system/include/visit_model_bond.h
%include configuration/include/site.h
%include configuration/include/particle.h
%include configuration/include/model_params.h
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
%include system/include/potential.h
%include system/include/potential_factory.h
%include system/include/system.h
%include system/include/utils.h
%include monte_carlo/include/trial_select.h
%include monte_carlo/include/trial_select_bond.h
%include monte_carlo/include/trial_select_angle.h
%include monte_carlo/include/perturb.h
%include monte_carlo/include/perturb_move.h
%include monte_carlo/include/perturb_translate.h
%include monte_carlo/include/perturb_distance.h
%include cluster/include/perturb_point_reflect.h
%include monte_carlo/include/perturb_volume.h
%include monte_carlo/include/trial_select_particle.h
%include cluster/include/select_particle_avb.h
%include cluster/include/select_particle_avb_divalent.h
%include cluster/include/select_cluster.h
%include monte_carlo/include/criteria.h
%include monte_carlo/include/metropolis.h
%include monte_carlo/include/trial_stage.h
%include monte_carlo/include/trial_compute.h
%include monte_carlo/include/trial.h
%include cluster/include/trial_avb4.h
%include monte_carlo/include/trial_factory.h
%include cluster/include/trial_avb2.h
%include cluster/include/trial_transfer_avb_divalent.h
%include cluster/include/trial_transfer_avb.h
%include monte_carlo/include/analyze.h
%include monte_carlo/include/analyze_factory.h
%include steppers/include/log_and_movie.h
%include steppers/include/seek_analyze.h
%include steppers/include/wall_clock_limit.h
%include steppers/include/profile_trials.h
%include steppers/include/volume.h
%include steppers/include/check_physicality.h
%include steppers/include/check.h
%include steppers/include/extensive_moments.h
%include steppers/include/cpu_time.h
%include steppers/include/mean_squared_displacement.h
%include steppers/include/num_particles.h
%include steppers/include/log.h
%include steppers/include/energy.h
%include steppers/include/criteria_writer.h
%include monte_carlo/include/modify.h
%include steppers/include/pair_distribution.h
%include steppers/include/check_energy.h
%include monte_carlo/include/modify_factory.h
%include monte_carlo/include/monte_carlo.h
%include steppers/include/check_energy_and_tune.h
%include steppers/include/seek_modify.h
%include steppers/include/tuner.h
%include steppers/include/check_properties.h
%include steppers/include/criteria_updater.h
%include monte_carlo/include/seek_num_particles.h
%include monte_carlo/include/trials.h
%include cluster/include/trial_rigid_cluster.h
%include monte_carlo/include/trial_compute_volume.h
%include monte_carlo/include/trial_compute_remove.h
%include cluster/include/compute_move_cluster.h
%include cluster/include/compute_add_avb.h
%include cluster/include/compute_remove_avb.h
%include cluster/include/compute_gca.h
%include monte_carlo/include/trial_compute_add.h
%include cluster/include/compute_add_avb_divalent.h
%include monte_carlo/include/trial_compute_move.h
%include monte_carlo/include/trial_move.h
%include cluster/include/compute_remove_avb_divalent.h
%include cluster/include/compute_avb2.h
%include cluster/include/compute_avb4.h
%include monte_carlo/include/constraint.h
%include monte_carlo/include/constrain_num_particles.h
%include configuration/include/visit_configuration.h
%include configuration/include/file_xyz.h
%include steppers/include/movie.h
%include configuration/include/utils.h
%include configuration/include/file_lmp.h
%include configuration/include/domain.h
%include math/include/random.h
%include math/include/constants.h
%include math/include/minimize.h
%include math/include/golden_search.h
%include math/include/quadratic_equation.h
%include math/include/formula_exponential.h
%include math/include/matrix.h
%include monte_carlo/include/perturb_rotate.h
%include cluster/include/perturb_rotate_com.h
%include cluster/include/perturb_move_avb.h
%include cluster/include/perturb_add_avb.h
%include monte_carlo/include/perturb_anywhere.h
%include monte_carlo/include/perturb_add.h
%include monte_carlo/include/perturb_remove.h
%include monte_carlo/include/perturb_distance_angle.h
%include math/include/random_mt19937.h
%include math/include/solver.h
%include math/include/solver_newton_raphson.h
%include math/include/solver_bisection.h
%include math/include/solver_brent_dekker.h
