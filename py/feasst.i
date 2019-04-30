/* This is an interface file for swig.
   This file is created by dev/tools/depend.py . Modifications to this
   file will likely be overwritten in the future. Thus, edit depend.py
   instead.

   usage: /path/to/feasst/py/run.sh /path/to/feasst/dev/tools/depend.py -s /path/to/feasst
 */

%module feasst

%{
#include "core/include/file.h"
#include "core/include/accumulator.h"
#include "core/include/properties.h"
#include "core/include/typed_entity.h"
#include "core/include/arguments.h"
#include "core/include/position.h"
#include "core/include/tunable.h"
#include "core/include/formula.h"
#include "core/include/formula_exponential.h"
#include "core/include/histogram.h"
#include "core/include/formula_polynomial.h"
#include "core/include/custom_exception.h"
#include "core/include/utils.h"
#include "core/include/debug.h"
#include "core/include/constants.h"
#include "core/include/random.h"
#include "core/include/utils_math.h"
#include "core/include/physical_constants.h"
#include "core/include/bond.h"
#include "core/include/utils_io.h"
#include "core/include/site.h"
#include "core/include/particle.h"
#include "core/include/file_lmp.h"
#include "core/include/model_params.h"
#include "core/include/group.h"
#include "core/include/particle_factory.h"
#include "core/include/select.h"
#include "core/include/select_position.h"
#include "core/include/cells.h"
#include "core/include/domain.h"
#include "core/include/visit_particles.h"
#include "core/include/stepper.h"
#include "core/include/matrix.h"
#include "core/include/configuration.h"
#include "core/include/file_xyz.h"
#include "core/include/select_list.h"
#include "core/include/model.h"
#include "core/include/visit_model.h"
#include "core/include/visit_model_cell.h"
#include "core/include/long_range_corrections.h"
#include "core/include/model_three_body.h"
#include "patch/include/visit_model_inner_patch.h"
#include "core/include/model_two_body.h"
#include "models/include/model_yukawa.h"
#include "core/include/model_lj.h"
#include "models/include/model_lj_alpha.h"
#include "models/include/model_lj_cut_shift.h"
#include "models/include/model_lj_force_shift.h"
#include "core/include/model_square_well.h"
#include "core/include/model_hard_sphere.h"
#include "core/include/visit_model_intra.h"
#include "core/include/bond_visitor.h"
#include "core/include/model_one_body.h"
#include "core/include/model_empty.h"
#include "core/include/potential.h"
#include "core/include/potential_factory.h"
#include "core/include/system.h"
#include "core/include/criteria.h"
#include "core/include/criteria_mayer.h"
#include "core/include/criteria_metropolis.h"
#include "core/include/perturb.h"
#include "core/include/perturb_configs.h"
#include "core/include/perturb_move.h"
#include "core/include/perturb_rotate.h"
#include "core/include/perturb_transfer.h"
#include "core/include/trial.h"
#include "core/include/trial_move.h"
#include "core/include/trial_rotate.h"
#include "core/include/trial_transfer.h"
#include "core/include/rosenbluth.h"
#include "core/include/trial_factory.h"
#include "core/include/modify.h"
#include "core/include/check.h"
#include "core/include/wall_clock_limit.h"
#include "core/include/tuner.h"
#include "core/include/modify_factory.h"
#include "core/include/analyze.h"
#include "core/include/movie.h"
#include "core/include/analyze_factory.h"
#include "core/include/criteria_writer.h"
#include "core/include/log.h"
#include "core/include/monte_carlo.h"
#include "core/include/perturb_translate.h"
#include "core/include/trial_translate.h"
#include "core/include/visit_configuration.h"
#include "ewald/include/model_charge_intra.h"
#include "ewald/include/ewald.h"
#include "ewald/include/model_charge_self.h"
#include "ewald/include/model_charge_screened.h"
#include "flat_histogram/include/ln_probability_distribution.h"
#include "flat_histogram/include/bias.h"
#include "flat_histogram/include/bias_wang_landau.h"
#include "flat_histogram/include/macrostate.h"
#include "flat_histogram/include/collection_matrix.h"
#include "flat_histogram/include/bias_transition_matrix.h"
#include "flat_histogram/include/macrostate_num_particles.h"
#include "flat_histogram/include/macrostate_accumulator.h"
#include "flat_histogram/include/criteria_flat_histogram.h"
#include "confinement/include/shape.h"
#include "confinement/include/model_hard_shape.h"
#include "confinement/include/half_space.h"
#include "chain/include/trial_crankshaft.h"
#include "chain/include/perturb_regrow.h"
#include "chain/include/analyze_rigid_bonds.h"
#include "chain/include/trial_regrow.h"
#include "chain/include/trial_reptate.h"
#include "chain/include/trial_pivot.h"
#include "example/include/model_example.h"
using namespace feasst;
%}

%include "std_string.i"
%include "std_vector.i"
%include "std_shared_ptr.i"
%include "std_iostream.i"
%template(IntVector) std::vector<int>;
%template(DoubleVector) std::vector<double>;
using namespace std;
%pythonnondynamic;
%include "std_map.i"
%template(args) std::map<std::string, std::string>;

%shared_ptr(feasst::Accumulator);
%shared_ptr(feasst::Properties);
%shared_ptr(feasst::PropertiedEntity);
%shared_ptr(feasst::TypedEntity);
%shared_ptr(feasst::Arguments);
%shared_ptr(feasst::Position);
%shared_ptr(feasst::PositionSpherical);
%shared_ptr(feasst::SpatialEntity);
%shared_ptr(feasst::Tunable);
%shared_ptr(feasst::Formula);
%shared_ptr(feasst::FormulaExponential);
%shared_ptr(feasst::Histogram);
%shared_ptr(feasst::FormulaPolynomial);
%shared_ptr(feasst::CustomException);
%shared_ptr(feasst::Random);
%shared_ptr(feasst::Bond);
%shared_ptr(feasst::Angle);
%shared_ptr(feasst::Dihedral);
%shared_ptr(feasst::Improper);
%shared_ptr(feasst::Site);
%shared_ptr(feasst::Particle);
%shared_ptr(feasst::FileLMP);
%shared_ptr(feasst::ModelParam);
%shared_ptr(feasst::Epsilon);
%shared_ptr(feasst::Sigma);
%shared_ptr(feasst::CutOff);
%shared_ptr(feasst::Charge);
%shared_ptr(feasst::ModelParams);
%shared_ptr(feasst::Group);
%shared_ptr(feasst::ParticleFactory);
%shared_ptr(feasst::Select);
%shared_ptr(feasst::SelectGroup);
%shared_ptr(feasst::SelectPosition);
%shared_ptr(feasst::Cells);
%shared_ptr(feasst::Domain);
%shared_ptr(feasst::VisitParticles);
%shared_ptr(feasst::LoopOneBody);
%shared_ptr(feasst::Stepper);
%shared_ptr(feasst::Matrix);
%shared_ptr(feasst::MatrixThreeByThree);
%shared_ptr(feasst::RotationMatrix);
%shared_ptr(feasst::Configuration);
%shared_ptr(feasst::FileVMD);
%shared_ptr(feasst::FileXYZ);
%shared_ptr(feasst::SelectList);
%shared_ptr(feasst::Model);
%shared_ptr(feasst::VisitModelInner);
%shared_ptr(feasst::VisitModel);
%shared_ptr(feasst::VisitModelCell);
%shared_ptr(feasst::LongRangeCorrections);
%shared_ptr(feasst::ModelThreeBody);
%shared_ptr(feasst::PatchAngle);
%shared_ptr(feasst::CosPatchAngle);
%shared_ptr(feasst::VisitModelInnerPatch);
%shared_ptr(feasst::ModelTwoBody);
%shared_ptr(feasst::ModelTwoBodyFactory);
%shared_ptr(feasst::ModelYukawa);
%shared_ptr(feasst::ModelLJ);
%shared_ptr(feasst::ModelLJAlpha);
%shared_ptr(feasst::EnergyAtCutoff);
%shared_ptr(feasst::ModelLJCutShift);
%shared_ptr(feasst::EnergyDerivAtCutoff);
%shared_ptr(feasst::ModelLJForceShift);
%shared_ptr(feasst::ModelSquareWell);
%shared_ptr(feasst::ModelHardSphere);
%shared_ptr(feasst::VisitModelIntra);
%shared_ptr(feasst::BondTwoBody);
%shared_ptr(feasst::BondSquareWell);
%shared_ptr(feasst::BondThreeBody);
%shared_ptr(feasst::AngleSquareWell);
%shared_ptr(feasst::BondVisitor);
%shared_ptr(feasst::ModelOneBody);
%shared_ptr(feasst::ModelEmpty);
%shared_ptr(feasst::Potential);
%shared_ptr(feasst::PotentialFactory);
%shared_ptr(feasst::System);
%shared_ptr(feasst::Criteria);
%shared_ptr(feasst::CriteriaMayer);
%shared_ptr(feasst::CriteriaMetropolis);
%shared_ptr(feasst::Perturb);
%shared_ptr(feasst::PerturbOptRevert);
%shared_ptr(feasst::PerturbConfigs);
%shared_ptr(feasst::PerturbSelectMove);
%shared_ptr(feasst::PerturbRotate);
%shared_ptr(feasst::PerturbAdd);
%shared_ptr(feasst::PerturbRemove);
%shared_ptr(feasst::Trial);
%shared_ptr(feasst::TrialMove);
%shared_ptr(feasst::TrialRotate);
%shared_ptr(feasst::TrialTransfer);
%shared_ptr(feasst::Rosenbluth);
%shared_ptr(feasst::Stage);
%shared_ptr(feasst::StageFactory);
%shared_ptr(feasst::TrialFactory);
%shared_ptr(feasst::Modify);
%shared_ptr(feasst::ModifyUpdateOnly);
%shared_ptr(feasst::Check);
%shared_ptr(feasst::EnergyCheck);
%shared_ptr(feasst::WallClockLimit);
%shared_ptr(feasst::Tuner);
%shared_ptr(feasst::ModifyFactory);
%shared_ptr(feasst::Analyze);
%shared_ptr(feasst::AnalyzeWriteOnly);
%shared_ptr(feasst::AnalyzeUpdateOnly);
%shared_ptr(feasst::Movie);
%shared_ptr(feasst::AnalyzeFactory);
%shared_ptr(feasst::CriteriaWriter);
%shared_ptr(feasst::Log);
%shared_ptr(feasst::MonteCarlo);
%shared_ptr(feasst::PerturbTranslate);
%shared_ptr(feasst::TrialTranslate);
%shared_ptr(feasst::StagedTrial);
%shared_ptr(feasst::TrialStagedTranslate);
%shared_ptr(feasst::VisitConfiguration);
%shared_ptr(feasst::ModelChargeIntra);
%shared_ptr(feasst::Ewald);
%shared_ptr(feasst::ModelChargeSelf);
%shared_ptr(feasst::ModelChargeScreened);
%shared_ptr(feasst::LnProbabilityDistribution);
%shared_ptr(feasst::Bias);
%shared_ptr(feasst::BiasWangLandau);
%shared_ptr(feasst::Macrostate);
%shared_ptr(feasst::TripleBandedCollectionMatrix);
%shared_ptr(feasst::BiasTransitionMatrix);
%shared_ptr(feasst::MacrostateNumParticles);
%shared_ptr(feasst::MacrostateAccumulator);
%shared_ptr(feasst::MacrostateAccumulatorFactory);
%shared_ptr(feasst::BinEnergy);
%shared_ptr(feasst::CriteriaFlatHistogram);
%shared_ptr(feasst::Shape);
%shared_ptr(feasst::ShapedEntity);
%shared_ptr(feasst::ShapeIntersect);
%shared_ptr(feasst::ModelHardShape);
%shared_ptr(feasst::HalfSpace);
%shared_ptr(feasst::TrialCrankshaft);
%shared_ptr(feasst::PerturbRegrow);
%shared_ptr(feasst::AnalyzeRigidBonds);
%shared_ptr(feasst::StageFactoryRegrow);
%shared_ptr(feasst::TrialRegrow);
%shared_ptr(feasst::TrialReptate);
%shared_ptr(feasst::TrialPivot);
%shared_ptr(feasst::ModelExample);
%include core/include/file.h
%include core/include/accumulator.h
%include core/include/properties.h
%include core/include/typed_entity.h
%include core/include/arguments.h
%include core/include/position.h
%include core/include/tunable.h
%include core/include/formula.h
%include core/include/formula_exponential.h
%include core/include/histogram.h
%include core/include/formula_polynomial.h
%include core/include/custom_exception.h
%include core/include/utils.h
%include core/include/debug.h
%include core/include/constants.h
%include core/include/random.h
%include core/include/utils_math.h
%include core/include/physical_constants.h
%include core/include/bond.h
%include core/include/utils_io.h
%include core/include/site.h
%include core/include/particle.h
%include core/include/file_lmp.h
%include core/include/model_params.h
%include core/include/group.h
%include core/include/particle_factory.h
%include core/include/select.h
%include core/include/select_position.h
%include core/include/cells.h
%include core/include/domain.h
%include core/include/visit_particles.h
%include core/include/stepper.h
%include core/include/matrix.h
%include core/include/configuration.h
%include core/include/file_xyz.h
%include core/include/select_list.h
%include core/include/model.h
%include core/include/visit_model.h
%include core/include/visit_model_cell.h
%include core/include/long_range_corrections.h
%include core/include/model_three_body.h
%include patch/include/visit_model_inner_patch.h
%include core/include/model_two_body.h
%include models/include/model_yukawa.h
%include core/include/model_lj.h
%include models/include/model_lj_alpha.h
%include models/include/model_lj_cut_shift.h
%include models/include/model_lj_force_shift.h
%include core/include/model_square_well.h
%include core/include/model_hard_sphere.h
%include core/include/visit_model_intra.h
%include core/include/bond_visitor.h
%include core/include/model_one_body.h
%include core/include/model_empty.h
%include core/include/potential.h
%include core/include/potential_factory.h
%include core/include/system.h
%include core/include/criteria.h
%include core/include/criteria_mayer.h
%include core/include/criteria_metropolis.h
%include core/include/perturb.h
%include core/include/perturb_configs.h
%include core/include/perturb_move.h
%include core/include/perturb_rotate.h
%include core/include/perturb_transfer.h
%include core/include/trial.h
%include core/include/trial_move.h
%include core/include/trial_rotate.h
%include core/include/trial_transfer.h
%include core/include/rosenbluth.h
%include core/include/trial_factory.h
%include core/include/modify.h
%include core/include/check.h
%include core/include/wall_clock_limit.h
%include core/include/tuner.h
%include core/include/modify_factory.h
%include core/include/analyze.h
%include core/include/movie.h
%include core/include/analyze_factory.h
%include core/include/criteria_writer.h
%include core/include/log.h
%include core/include/monte_carlo.h
%include core/include/perturb_translate.h
%include core/include/trial_translate.h
%include core/include/visit_configuration.h
%include ewald/include/model_charge_intra.h
%include ewald/include/ewald.h
%include ewald/include/model_charge_self.h
%include ewald/include/model_charge_screened.h
%include flat_histogram/include/ln_probability_distribution.h
%include flat_histogram/include/bias.h
%include flat_histogram/include/bias_wang_landau.h
%include flat_histogram/include/macrostate.h
%include flat_histogram/include/collection_matrix.h
%include flat_histogram/include/bias_transition_matrix.h
%include flat_histogram/include/macrostate_num_particles.h
%include flat_histogram/include/macrostate_accumulator.h
%include flat_histogram/include/criteria_flat_histogram.h
%include confinement/include/shape.h
%include confinement/include/model_hard_shape.h
%include confinement/include/half_space.h
%include chain/include/trial_crankshaft.h
%include chain/include/perturb_regrow.h
%include chain/include/analyze_rigid_bonds.h
%include chain/include/trial_regrow.h
%include chain/include/trial_reptate.h
%include chain/include/trial_pivot.h
%include example/include/model_example.h
