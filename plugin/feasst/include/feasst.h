#include "math/include/constants.h"
#include "math/include/euler.h"
#include "math/include/formula.h"
#include "math/include/formula_exponential.h"
#include "math/include/formula_polynomial.h"
#include "math/include/histogram.h"
#include "math/include/matrix.h"
#include "math/include/minimize.h"
#include "math/include/golden_search.h"
#include "math/include/position.h"
#include "math/include/quadratic_equation.h"
#include "math/include/random.h"
#include "math/include/random_modulo.h"
#include "math/include/random_mt19937.h"
#include "math/include/solver.h"
#include "math/include/solver_bisection.h"
#include "math/include/solver_brent_dekker.h"
#include "math/include/solver_newton_raphson.h"
#include "math/include/table.h"
#include "math/include/utils_math.h"
#include "math/include/accumulator.h"
#include "monte_carlo/include/acceptance.h"
#include "monte_carlo/include/action.h"
#include "example/include/action_example.h"
#include "monte_carlo/include/constraint.h"
#include "monte_carlo/include/constrain_volume_by_cutoff.h"
#include "monte_carlo/include/constrain_num_particles.h"
#include "monte_carlo/include/perturb.h"
#include "monte_carlo/include/perturb_move.h"
#include "monte_carlo/include/perturb_remove.h"
#include "monte_carlo/include/perturb_rotate.h"
#include "monte_carlo/include/perturb_translate.h"
#include "monte_carlo/include/perturb_anywhere.h"
#include "monte_carlo/include/perturb_add.h"
#include "monte_carlo/include/perturb_volume.h"
#include "monte_carlo/include/ref_potential.h"
#include "monte_carlo/include/remove_analyze.h"
#include "monte_carlo/include/remove_modify.h"
#include "monte_carlo/include/remove_trial.h"
#include "monte_carlo/include/trial_compute.h"
#include "monte_carlo/include/trial_compute_add.h"
#include "monte_carlo/include/trial_compute_move.h"
#include "monte_carlo/include/trial_compute_remove.h"
#include "monte_carlo/include/trial_compute_translate.h"
#include "monte_carlo/include/trial_compute_volume.h"
#include "monte_carlo/include/trial_select.h"
#include "monte_carlo/include/trial_select_all.h"
#include "monte_carlo/include/trial_select_bond.h"
#include "monte_carlo/include/trial_select_angle.h"
#include "monte_carlo/include/trial_select_dihedral.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/tunable.h"
#include "monte_carlo/include/write_checkpoint.h"
#include "monte_carlo/include/write_model_params.h"
#include "monte_carlo/include/convert_to_ref_potential.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/run.h"
#include "monte_carlo/include/stepper.h"
#include "monte_carlo/include/analyze.h"
#include "monte_carlo/include/analyze_factory.h"
#include "example/include/analyze_example.h"
#include "monte_carlo/include/modify.h"
#include "monte_carlo/include/modify_factory.h"
#include "monte_carlo/include/remove.h"
#include "chain/include/perturb_connector.h"
#include "chain/include/perturb_crankshaft.h"
#include "chain/include/perturb_crankshaft_small.h"
#include "chain/include/perturb_library.h"
#include "chain/include/perturb_particle_pivot.h"
#include "chain/include/perturb_pivot.h"
#include "chain/include/perturb_position_swap.h"
#include "chain/include/perturb_site_type.h"
#include "chain/include/perturb_to_anchor.h"
#include "chain/include/select_branch.h"
#include "chain/include/select_crankshaft_small.h"
#include "chain/include/select_particle_pivot.h"
#include "chain/include/select_perturbed.h"
#include "chain/include/select_segment.h"
#include "chain/include/select_end_segment.h"
#include "chain/include/select_site_of_type.h"
#include "chain/include/select_two_sites.h"
#include "chain/include/end_to_end_distance.h"
#include "chain/include/radius_of_gyration.h"
#include "morph/include/perturb_particle_type.h"
#include "prefetch/include/prefetch.h"
#include "server/include/server.h"
#include "server/include/listen.h"
#include "steppers/include/seek_analyze.h"
#include "steppers/include/seek_modify.h"
#include "steppers/include/write_stepper.h"
#include "steppers/include/check.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/check_physicality.h"
#include "steppers/include/check_properties.h"
#include "steppers/include/chirality_2d.h"
#include "steppers/include/criteria_updater.h"
#include "steppers/include/criteria_writer.h"
#include "steppers/include/density.h"
#include "steppers/include/density_profile.h"
#include "steppers/include/energy.h"
#include "steppers/include/extensive_moments.h"
#include "steppers/include/ghost_trial_volume.h"
#include "steppers/include/heat_capacity.h"
#include "steppers/include/increment_phase.h"
#include "steppers/include/movie.h"
#include "steppers/include/num_particles.h"
#include "steppers/include/profile_cpu.h"
#include "steppers/include/profile_trials.h"
#include "steppers/include/scattering.h"
#include "steppers/include/tune.h"
#include "steppers/include/volume.h"
#include "steppers/include/wall_clock_limit.h"
#include "steppers/include/wrap_particles.h"
#include "system/include/bond_four_body.h"
#include "models/include/dihedral_trappe.h"
#include "models/include/dihedral_ryckaert_bellemans.h"
#include "models/include/dihedral_harmonic.h"
#include "system/include/bond_three_body.h"
#include "system/include/angle_square_well.h"
#include "models/include/angle_harmonic.h"
#include "system/include/bond_two_body.h"
#include "models/include/fene.h"
#include "system/include/bond_square_well.h"
#include "models/include/bond_harmonic.h"
#include "system/include/bond_visitor.h"
#include "steppers/include/log.h"
#include "system/include/model.h"
#include "system/include/model_one_body.h"
#include "system/include/model_empty.h"
#include "system/include/model_three_body.h"
#include "system/include/model_two_body.h"
#include "system/include/lennard_jones.h"
#include "server/include/model_server.h"
#include "models/include/f3c.h"
#include "models/include/yukawa.h"
#include "models/include/two_body_table.h"
#include "models/include/table_potential.h"
#include "models/include/square_well.h"
#include "system/include/ideal_gas.h"
#include "system/include/hard_sphere.h"
#include "example/include/model_example.h"
#include "system/include/model_two_body_factory.h"
#include "system/include/model_two_body_table.h"
#include "system/include/potential_factory.h"
#include "system/include/rigid_angle.h"
#include "system/include/rigid_bond.h"
#include "monte_carlo/include/perturb_distance.h"
#include "chain/include/perturb_reptate.h"
#include "monte_carlo/include/perturb_distance_angle.h"
#include "chain/include/perturb_distance_angle_connector.h"
#include "chain/include/perturb_branch.h"
#include "system/include/rigid_dihedral.h"
#include "chain/include/analyze_bonds.h"
#include "monte_carlo/include/perturb_dihedral.h"
#include "system/include/synchronize_data.h"
#include "system/include/energy_map.h"
#include "monte_carlo/include/criteria.h"
#include "mayer/include/mayer_sampling.h"
#include "monte_carlo/include/always_reject.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial.h"
#include "morph/include/trial_morph_expanded.h"
#include "morph/include/trial_morph.h"
#include "chain/include/ghost_trial_grow.h"
#include "chain/include/trial_swap_sites.h"
#include "chain/include/trial_grow_linear.h"
#include "monte_carlo/include/trial_factory.h"
#include "chain/include/trial_grow.h"
#include "monte_carlo/include/trial_transfer.h"
#include "monte_carlo/include/trial_volume.h"
#include "monte_carlo/include/trial_remove.h"
#include "monte_carlo/include/trial_move.h"
#include "chain/include/trial_reptate.h"
#include "chain/include/trial_pivot.h"
#include "chain/include/trial_particle_pivot.h"
#include "chain/include/trial_crankshaft_small.h"
#include "chain/include/trial_crankshaft.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/trial_rotate.h"
#include "monte_carlo/include/trial_add.h"
#include "system/include/thermo_params.h"
#include "system/include/visit_model.h"
#include "system/include/long_range_corrections.h"
#include "opt_lj/include/visit_model_opt_rpm.h"
#include "opt_lj/include/visit_model_opt_lj.h"
#include "system/include/dont_visit_model.h"
#include "system/include/visit_model_bond.h"
#include "system/include/visit_model_cell.h"
#include "system/include/visit_model_inner.h"
#include "server/include/visit_model_inner_server.h"
#include "aniso/include/visit_model_inner_table.h"
#include "aniso/include/visit_model_inner_nn.h"
#include "system/include/visit_model_intra.h"
#include "system/include/visit_model_intra_map.h"
#include "system/include/potential.h"
#include "system/include/system.h"
#include "morph/include/compute_morph.h"
#include "system/include/visit_model_cutoff_outer.h"
#include "threads/include/thread.h"
#include "threads/include/thread_omp.h"
#include "utils/include/argument_parse.h"
#include "utils/include/arguments.h"
#include "utils/include/arguments_extra.h"
#include "utils/include/cache.h"
#include "utils/include/custom_exception.h"
#include "utils/include/debug.h"
#include "utils/include/file.h"
#include "utils/include/max_precision.h"
#include "utils/include/io.h"
#include "utils/include/progress_report.h"
#include "utils/include/timer.h"
#include "steppers/include/cpu_time.h"
#include "utils/include/utils.h"
#include "utils/include/checkpoint.h"
#include "utils/include/serialize.h"
#include "utils/include/serialize_extra.h"
#include "utils/include/timer_rdtsc.h"
#include "gibbs/include/trial_gibbs_particle_transfer.h"
#include "gibbs/include/check_constant_volume.h"
#include "gibbs/include/compute_gibbs_particle_transfer.h"
#include "gibbs/include/compute_gibbs_volume_transfer.h"
#include "gibbs/include/copy_following_lines.h"
#include "gibbs/include/copy_next_line.h"
#include "gibbs/include/end_copy.h"
#include "gibbs/include/gibbs_initialize.h"
#include "gibbs/include/trial_gibbs_volume_transfer.h"
#include "flat_histogram/include/bias.h"
#include "flat_histogram/include/collection_matrix.h"
#include "flat_histogram/include/flat_histogram.h"
#include "flat_histogram/include/ln_probability.h"
#include "flat_histogram/include/macrostate.h"
#include "morph/include/macrostate_morph.h"
#include "flat_histogram/include/ensemble.h"
#include "flat_histogram/include/macrostate_energy.h"
#include "flat_histogram/include/macrostate_num_particles.h"
#include "flat_histogram/include/macrostate_position.h"
#include "flat_histogram/include/transition_matrix.h"
#include "flat_histogram/include/wang_landau.h"
#include "flat_histogram/include/window.h"
#include "flat_histogram/include/collection_matrix_splice.h"
#include "flat_histogram/include/clones.h"
#include "flat_histogram/include/window_custom.h"
#include "flat_histogram/include/window_exponential.h"
#include "flat_histogram/include/wltm.h"
#include "beta_expanded/include/compute_beta.h"
#include "beta_expanded/include/macrostate_beta.h"
#include "beta_expanded/include/perturb_beta.h"
#include "beta_expanded/include/select_nothing.h"
#include "beta_expanded/include/trial_beta.h"
#include "charge/include/charge_screened.h"
#include "charge/include/charge_screened_intra.h"
#include "charge/include/charge_self.h"
#include "charge/include/compute_add_multiple.h"
#include "charge/include/compute_remove_multiple.h"
#include "charge/include/debye_huckel.h"
#include "charge/include/electric_field.h"
#include "charge/include/slab_correction.h"
#include "charge/include/trial_add_multiple.h"
#include "charge/include/trial_remove_multiple.h"
#include "charge/include/trial_transfer_multiple.h"
#include "charge/include/utils.h"
#include "charge/include/coulomb.h"
#include "cluster/include/compute_add_avb.h"
#include "cluster/include/compute_avb2.h"
#include "cluster/include/compute_avb4.h"
#include "cluster/include/compute_gca.h"
#include "cluster/include/compute_move_cluster.h"
#include "cluster/include/compute_remove_avb.h"
#include "cluster/include/perturb_move_avb.h"
#include "cluster/include/perturb_add_avb.h"
#include "cluster/include/perturb_point_reflect.h"
#include "cluster/include/perturb_rotate_com.h"
#include "cluster/include/select_cluster.h"
#include "cluster/include/trial_add_avb.h"
#include "cluster/include/trial_add_avb_divalent.h"
#include "cluster/include/trial_avb2.h"
#include "cluster/include/trial_avb4.h"
#include "cluster/include/trial_remove_avb.h"
#include "cluster/include/trial_remove_avb_divalent.h"
#include "cluster/include/trial_rigid_cluster.h"
#include "cluster/include/trial_rotate_cluster.h"
#include "cluster/include/trial_transfer_avb.h"
#include "cluster/include/trial_transfer_avb_divalent.h"
#include "cluster/include/trial_translate_cluster.h"
#include "cluster/include/analyze_cluster.h"
#include "cluster/include/calculate_cluster.h"
#include "configuration/include/configuration.h"
#include "charge/include/ewald.h"
#include "charge/include/check_net_charge.h"
#include "configuration/include/domain.h"
#include "configuration/include/neighbor_criteria.h"
#include "cluster/include/energy_map_neighbor.h"
#include "cluster/include/energy_map_neighbor_criteria.h"
#include "cluster/include/energy_map_all.h"
#include "cluster/include/energy_map_all_criteria.h"
#include "configuration/include/physical_constants.h"
#include "configuration/include/properties.h"
#include "configuration/include/model_params.h"
#include "steppers/include/pair_distribution.h"
#include "aniso/include/anisotropic.h"
#include "example/include/model_param_example.h"
#include "system/include/cutoff_outer.h"
#include "models/include/jagla.h"
#include "models/include/lennard_jones_alpha.h"
#include "models/include/lennard_jones_cut_shift.h"
#include "models/include/lennard_jones_force_shift.h"
#include "models/include/mie_parameters.h"
#include "models/include/mie.h"
#include "models/include/two_body_alpha.h"
#include "configuration/include/group.h"
#include "configuration/include/bond.h"
#include "configuration/include/select.h"
#include "steppers/include/mean_squared_displacement.h"
#include "chain/include/select_reptate.h"
#include "monte_carlo/include/rosenbluth.h"
#include "cluster/include/select_particle_avb_divalent.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/compute_remove_avb_divalent.h"
#include "system/include/cells.h"
#include "cluster/include/compute_add_avb_divalent.h"
#include "configuration/include/site.h"
#include "configuration/include/particle.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/file_particle.h"
#include "configuration/include/visit_particles.h"
#include "configuration/include/visit_configuration.h"
#include "configuration/include/file_xyz.h"
#include "steppers/include/read_config_from_file.h"
#include "aniso/include/rotator.h"
#include "aniso/include/tabulate_two_rigid_body_3D.h"
#include "configuration/include/file_vmd.h"
#include "aniso/include/backmap.h"
#include "confinement/include/model_table_cartesian.h"
#include "confinement/include/trial_anywhere.h"
#include "confinement/include/background.h"
#include "confinement/include/henry_coefficient.h"
#include "confinement/include/zero_background.h"
#include "egce/include/a_equal_b.h"
#include "egce/include/a_half_b.h"
#include "model_expanded/include/compute_model.h"
#include "model_expanded/include/constrain_model_index.h"
#include "model_expanded/include/macrostate_model.h"
#include "model_expanded/include/model_expanded.h"
#include "model_expanded/include/perturb_model.h"
#include "model_expanded/include/trial_model.h"
#include "patch/include/patch_angle.h"
#include "patch/include/file_xyz_patch.h"
#include "patch/include/file_xyz_spherocylinder.h"
#include "patch/include/movie_patch.h"
#include "patch/include/movie_spherocylinder.h"
#include "patch/include/solid_of_revolution_table.h"
#include "patch/include/spherocylinder.h"
#include "patch/include/two_particle_contact.h"
#include "patch/include/visit_model_inner_patch.h"
#include "shape/include/formula_sine_wave.h"
#include "shape/include/shape.h"
#include "shape/include/half_space_tilted.h"
#include "shape/include/half_space.h"
#include "shape/include/half_space_sine.h"
#include "shape/include/cylinder.h"
#include "shape/include/cuboid.h"
#include "confinement/include/model_square_well_shape.h"
#include "confinement/include/model_lj_shape.h"
#include "confinement/include/model_hard_shape.h"
#include "shape/include/shape_file.h"
#include "shape/include/shape_intersect.h"
#include "shape/include/finite_cylinder.h"
#include "shape/include/shape_union.h"
#include "shape/include/slab.h"
#include "shape/include/slab_sine.h"
#include "shape/include/sphere.h"
#include "shape/include/supertoroid.h"
std::shared_ptr<feasst::ComputeBeta> __feasst__ComputeBeta = std::make_shared<feasst::ComputeBeta>();
std::shared_ptr<feasst::TrialGrow> __feasst__TrialGrow = std::make_shared<feasst::TrialGrow>();
std::shared_ptr<feasst::TrialParticlePivot> __feasst__TrialParticlePivot = std::make_shared<feasst::TrialParticlePivot>();
std::shared_ptr<feasst::ModelHardShape> __feasst__ModelHardShape = std::make_shared<feasst::ModelHardShape>();
std::shared_ptr<feasst::Sphere> __feasst__Sphere = std::make_shared<feasst::Sphere>();
std::shared_ptr<feasst::AEqualB> __feasst__AEqualB = std::make_shared<feasst::AEqualB>();
std::shared_ptr<feasst::Ewald> __feasst__Ewald = std::make_shared<feasst::Ewald>();
std::shared_ptr<feasst::ModelExample> __feasst__ModelExample = std::make_shared<feasst::ModelExample>();
std::shared_ptr<feasst::MayerSampling> __feasst__MayerSampling = std::make_shared<feasst::MayerSampling>();
std::shared_ptr<feasst::LennardJonesAlpha> __feasst__LennardJonesAlpha = std::make_shared<feasst::LennardJonesAlpha>();
std::shared_ptr<feasst::ComputeMorph> __feasst__ComputeMorph = std::make_shared<feasst::ComputeMorph>();
std::shared_ptr<feasst::VisitModelInnerPatch> __feasst__VisitModelInnerPatch = std::make_shared<feasst::VisitModelInnerPatch>();
std::shared_ptr<feasst::VisitModelInnerTable> __feasst__VisitModelInnerTable = std::make_shared<feasst::VisitModelInnerTable>();
std::shared_ptr<feasst::ComputeGibbsParticleTransfer> __feasst__ComputeGibbsParticleTransfer = std::make_shared<feasst::ComputeGibbsParticleTransfer>();
std::shared_ptr<feasst::MacrostateModel> __feasst__MacrostateModel = std::make_shared<feasst::MacrostateModel>();
std::shared_ptr<feasst::Listen> __feasst__Listen = std::make_shared<feasst::Listen>();
std::shared_ptr<feasst::Tune> __feasst__Tune = std::make_shared<feasst::Tune>();
