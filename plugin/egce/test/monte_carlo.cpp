#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "system/include/lennard_jones.h"
#include "system/include/long_range_corrections.h"
#include "system/include/hard_sphere.h"
#include "system/include/model_two_body_factory.h"
#include "system/include/dont_visit_model.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trials.h"
#include "monte_carlo/include/run.h"
#include "steppers/include/check_properties.h"
#include "steppers/include/cpu_time.h"
#include "steppers/include/criteria_updater.h"
#include "steppers/include/criteria_writer.h"
#include "steppers/include/energy.h"
#include "steppers/include/check_energy_and_tune.h"
#include "steppers/include/log_and_movie.h"
#include "flat_histogram/include/flat_histogram.h"
#include "flat_histogram/include/macrostate_num_particles.h"
#include "flat_histogram/include/transition_matrix.h"
#include "flat_histogram/include/wang_landau.h"
#include "charge/include/check_net_charge.h"
#include "charge/include/utils.h"
#include "charge/include/charge_screened.h"
#include "cluster/include/energy_map_all.h"
#include "cluster/include/energy_map_neighbor.h"
#include "cluster/include/energy_map_neighbor_criteria.h"
#include "cluster/include/trial_transfer_avb.h"
#include "cluster/include/trial_transfer_avb_divalent.h"
#include "chain/include/trial_grow.h"
#include "morph/include/macrostate_morph.h"
#include "morph/include/trial_morph_expanded.h"
#include "egce/include/a_equal_b.h"
#include "egce/include/a_half_b.h"
#include "opt_lj/include/visit_model_opt_rpm.h"

namespace feasst {

MonteCarlo rpm_egce(const int min = 0,
  const std::string dual_cut = "-1",
  const int steps_per = 1e5,
  const bool energy = true,
  const int max = 4) {
  MonteCarlo mc;
  // mc.set(MakeRandomMT19937({{"seed", "default"}}));
  mc.set(rpm({
    {"cubic_box_length", "12"},
    {"cutoff", "4.891304347826090"},
    {"kmax_squared", "38"},
    {"alpha", str(6.87098396396261/12)},
    {"dual_cut", dual_cut}}));
  const double temperature = 0.047899460618081;
  const double beta_mu = -13.94;
  mc.set(MakeThermoParams({{"beta", str(1/temperature)},
    {"chemical_potential0", str(beta_mu*temperature)},
    {"chemical_potential1", str(beta_mu*temperature)}}));
  if (min == 1) {
    mc.get_system()->get_configuration()->add_particle_of_type(0);
  } else if (min > 1) {
    FATAL("not implemented for min > 1");
  }
  auto criteria = MakeFlatHistogram(
    MakeMacrostateNumParticles(
      Histogram({{"width", "1"}, {"max", str(max)}, {"min", str(min)}})),
    MakeTransitionMatrix({{"min_sweeps", "100"}}),
    MakeAEqualB({{"extra_A", "1"}}));
  mc.set(criteria);
  mc.add(MakeCriteriaUpdater({{"steps_per", str(steps_per)}}));
  mc.add(MakeCriteriaWriter({{"steps_per", str(steps_per)},
                             {"file_name", "tmp/rpm_egce_crit.txt"}}));
  mc.add(MakeLogAndMovie({{"steps_per", str(steps_per)}, {"file_name", "tmp/rpm_egce"}}));
  mc.add(MakeCheckEnergyAndTune({{"steps_per", str(steps_per)}, {"tolerance", "1e-8"}}));
  // mc.add(MakeCheckProperties({{"steps_per", str(steps_per)}}));
  // mc.add(MakeCPUTime({{"steps_per", str(5*steps_per)}}));
  mc.add(MakeCheckNetCharge({{"maximum", "1."}, {"minimum", str(-NEAR_ZERO)}}));
  if (energy) {
    mc.add(MakeEnergy({
      {"file_name", "tmp/rpm_egce_energy"},
      {"steps_per_update", "1"},
      {"steps_per_write", str(steps_per)},
      {"multistate", "true"}}));
  }
  return mc;
}

TEST(MonteCarlo, rpm_egce_fh_LONG) {
  for (int num_steps : {1, 2}) {
    INFO("num_steps: " << num_steps);
    MonteCarlo mc;
    std::string ref = "-1";
    if (num_steps == 1) {
      mc = rpm_egce();
    } else {
      mc = rpm_egce(0, "1.");
      ref = "0";
    }
    // mc.set(MakeRandomMT19937({{"seed", "1234"}}));
    // mc.set(MakeRandomMT19937({{"seed", "1603117667"}}));
    mc.add(MakeTrialTranslate({{"weight", "0.25"}, {"tunable_param", "0.1"}}));
    mc.add(MakeTrialTransfer({{"weight", "1."}, {"particle_type", "0"}}));
    mc.add(MakeTrialTransfer({
      {"weight", "1."},
      {"weight", "1."},
      {"particle_type", "1"},
      {"reference_index", ref},
      {"num_steps", str(num_steps)}}));
    MonteCarlo mc2 = test_serialize(mc);
    mc2.run_until_complete();

    std::stringstream ss;
    mc2.criteria().serialize(ss);
    FlatHistogram fh(ss);
    LnProbability lnpi3 = fh.bias().ln_prob().reduce(2);
    EXPECT_NEAR(lnpi3.value(0), -1.2994315780357, 0.11);
    EXPECT_NEAR(lnpi3.value(1), -1.08646312498868, 0.11);
    EXPECT_NEAR(lnpi3.value(2), -0.941850889679828, 0.11);
    const std::vector<std::shared_ptr<Analyze> >& en = mc2.analyzers().back()->analyzers();
    EXPECT_NEAR(en[0]->accumulator().average(), 0, 1e-14);
    EXPECT_NEAR(en[1]->accumulator().average(), -0.115474, 1e-6);
    EXPECT_NEAR(en[2]->accumulator().average(), -0.939408, 0.02);
    EXPECT_NEAR(en[3]->accumulator().average(), -1.32485, 0.03);
    EXPECT_NEAR(en[4]->accumulator().average(), -2.02625, 0.05);
  }
}

//TEST(MonteCarlo, rpm_egce_fh_liquid_LONG) {
//  for (int num_steps : {1}) {
//  //for (int num_steps : {4}) {
//  //for (int num_steps : {1, 2}) {
//    INFO("num_steps: " << num_steps);
//    MonteCarlo mc;
//    std::string ref = "-1";
//    if (num_steps == 1) {
//      mc = rpm_egce(100, "-1", 1e3, true, 105);
//    } else {
//      mc = rpm_egce(100, "1.", 1e3, true, 105);
//      ref = "0";
//    }
//    // mc.set(MakeRandomMT19937({{"seed", "1234"}}));
//    mc.add(MakeTrialTranslate({{"weight", "0.25"}, {"tunable_param", "0.1"}}));
//    mc.add(MakeTrialTransfer({{"weight", "1."}, {"particle_type", "0"}}));
//    mc.add(MakeTrialTransfer({
//      {"weight", "1."},
//      {"weight", "1."},
//      {"particle_type", "1"},
//      {"reference_index", ref},
//      {"num_steps", str(num_steps)}}));
////    SeeeekNumParticles(100).with_metropolis(
////      MakeAEqualB({{"extra_A", "1"}}),
////      {{"beta", "0.01"}, {"chemical_potential0", "1"}, {"chemical_potential1", "1"}}).run(&mc);
//    MonteCarlo mc2 = test_serialize(mc);
//    mc2.run_until_complete();
//  }
//}

TEST(MonteCarlo, rpm_egce_fh_min1_VERY_LONG) {
  MonteCarlo mc = rpm_egce(1);
  mc.add(MakeTrialTranslate({{"weight", "0.25"}, {"tunable_param", "0.1"}}));
  mc.add(MakeTrialTransfer({{"weight", "1."}, {"particle_type", "0"}}));
  mc.add(MakeTrialTransfer({{"weight", "1."}, {"particle_type", "1"}}));
  mc.run_until_complete();
  test_serialize(mc);

  std::stringstream ss;
  mc.criteria().serialize(ss);
  FlatHistogram fh(ss);
  const LnProbability& lnpi = fh.bias().ln_prob();
  EXPECT_NEAR(lnpi.value(0), -5.05, 0.12);
  EXPECT_NEAR(lnpi.value(1), -0.77327, 0.12);
  EXPECT_NEAR(lnpi.value(2), -3.55107, 0.12);
  EXPECT_NEAR(lnpi.value(3), -0.686417, 0.12);
  const std::vector<std::shared_ptr<Analyze> >& en = mc.analyzers().back()->analyzers();
  EXPECT_NEAR(en[0]->accumulator().average(), -0.115474, 1e-6);
  EXPECT_NEAR(en[1]->accumulator().average(), -0.939408, 0.03);
  EXPECT_NEAR(en[2]->accumulator().average(), -1.32485, 0.04);
  EXPECT_NEAR(en[3]->accumulator().average(), -2.02625, 0.05);
}

TEST(MonteCarlo, rpm_egce_avb_fh_LONG) {
  MonteCarlo mc = rpm_egce(1);
  mc.add(MakeTrialTranslate({{"weight", "0.25"}, {"tunable_param", "0.1"}}));
  mc.set(MakeRandomMT19937({{"seed", "123"}}));
  // mc.set(MakeRandomMT19937({{"seed", "1346867550"}}));
  mc.add(MakeNeighborCriteria({{"maximum_distance", "3"},
                               {"minimum_distance", "1"},
                               {"site_type0", "0"},
                               {"site_type1", "1"},
                               {"potential_index", "1"}}));
  mc.set(1, MakePotential(MakeModelTwoBodyFactory({MakeHardSphere(),
                                                   MakeChargeScreened()}),
                      //MakeVisitModelOptRPM(MakeVisitModelInner(MakeEnergyMapNeighborCriteria(neighbor_criteria)))));
                      //MakeVisitModelOptRPM(MakeVisitModelInner(MakeEnergyMapNeighbor()))));
                      //MakeVisitModelOptRPM(MakeVisitModelInner(MakeEnergyMapAll()))));
                      MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
  mc.get_system()->energy();
  mc.add(MakeTrialTransferAVB({
    {"neighbor_index", "0"},
    {"weight", "1."},
    {"particle_type", "0"},
    {"target_particle_type", "1"}}));
  mc.add(MakeTrialTransferAVB({
    {"neighbor_index", "0"},
    {"weight", "1."},
    {"particle_type", "1"},
    {"target_particle_type", "0"}}));
  mc.run_until_complete();
  test_serialize(mc);

  std::stringstream ss;
  mc.criteria().serialize(ss);
  FlatHistogram fh(ss);
  const LnProbability& lnpi = fh.bias().ln_prob();
  EXPECT_NEAR(lnpi.value(0), -5.05, 0.1);
  EXPECT_NEAR(lnpi.value(1), -0.77327, 0.1);
  EXPECT_NEAR(lnpi.value(2), -3.55107, 0.1);
  EXPECT_NEAR(lnpi.value(3), -0.686417, 0.1);
  const std::vector<std::shared_ptr<Analyze> >& en = mc.analyzers().back()->analyzers();
  EXPECT_NEAR(en[0]->accumulator().average(), -0.115474, 1e-6);
  EXPECT_NEAR(en[1]->accumulator().average(), -0.939408, 0.02);
  EXPECT_NEAR(en[2]->accumulator().average(), -1.32485, 0.03);
  EXPECT_NEAR(en[3]->accumulator().average(), -2.02625, 0.045);
}

// HWH is there an issue with changing the charge in how Fourier is computed?
//TEST(MonteCarlo, rpm_growth_expanded_LONG) {
//  const int min = 0, steps_per = 1e0;
//  MonteCarlo mc = rpm_egce(min, "-1", steps_per, false);
//  mc.set(MakeRandomMT19937({{"seed", "123"}}));
//
//  // add growth expanded particle types with half charge
//  { Configuration * config = mc.get_system()->get_configuration();
//    const double q_plus = config->model_params().charge().value(0);
//    const double q_minus = config->model_params().charge().value(1);
//    config->add_particle_type(install_dir() +
//      "/plugin/charge/forcefield/rpm_plus.fstprt", "0.5");
//    config->add_particle_type(install_dir() +
//      "/plugin/charge/forcefield/rpm_minus.fstprt", "0.5");
//    config->set_model_param("charge", 2, 0.5*q_plus);
//    config->set_model_param("charge", 3, 0.5*q_minus);
//    config->set_model_param("cutoff", 2, 6);
//    config->set_model_param("cutoff", 3, 6);
//  }
//  const std::vector<std::vector<int> > grow_sequence = {{2, 3}, {0, 1}};
//  mc.set(MakeFlatHistogram(
//    MakeMacrostateMorph(
//      grow_sequence,
//      Histogram({{"width", str(1./grow_sequence.size())},
//                 {"max", "4"}, {"min", "0"}})),
//    MakeWangLandau({{"min_flatness", "25"}}),
//    // MakeTransitionMatrix({{"min_sweeps", "10"}}),
//    { {"beta", str(mc.criteria().beta())},
//      {"chemical_potential0", str(mc.criteria().chemical_potential(0))},
//      {"chemical_potential1", str(mc.criteria().chemical_potential(1))},
//      {"chemical_potential2", str(mc.criteria().chemical_potential(0))},
//      {"chemical_potential3", str(mc.criteria().chemical_potential(1))}}));
//  mc.initialize_criteria();
//  mc.initialize_analyzers();
//  mc.add(MakeTrialTranslate({{"weight", "0.25"}, {"tunable_param", "0.1"}}));
//  mc.add(MakeTrialMorphExpanded(grow_sequence));
//  mc.add(MakeEnergy({
//    {"file_name", "tmp/dival_egce_energy"},
//    {"steps_per_update", "1"},
//    {"steps_per_write", str(int(1e5))},
//    {"multistate", "true"}}));
//  mc.run_until_complete();
//  const LnProbability lnpi = FlatHistogram(mc.criteria()).bias().ln_prob();
//  EXPECT_NEAR(lnpi.value(0), -5.41, 0.1);
//  EXPECT_NEAR(lnpi.value(2), -3.42548, 0.1);
//  EXPECT_NEAR(lnpi.value(4), -2.02966, 0.1);
//  EXPECT_NEAR(lnpi.value(6), -1.35573, 0.1);
//  EXPECT_NEAR(lnpi.value(8), -1.16195, 0.1);
//  EXPECT_NEAR(lnpi.value(10), -1.34209, 0.1);
//  const std::vector<std::shared_ptr<Analyze> >& en =
//    mc.analyzers().back()->analyzers();
//  EXPECT_NEAR(en[0]->accumulator().average(), 0, 1e-14);
//  EXPECT_NEAR(en[2]->accumulator().average(), -1.462473, 0.03);
//  EXPECT_NEAR(en[4]->accumulator().average(), -3.015735, 0.05);
//  EXPECT_NEAR(en[6]->accumulator().average(), -4.879190, 0.08);
//  EXPECT_NEAR(en[8]->accumulator().average(), -6.829874, 0.12);
//  EXPECT_NEAR(en[10]->accumulator().average(), -8.815852, 0.16);
//}

MonteCarlo dival_egce(
    const int min = 0,
    const int max = 15,
    const int steps_per = 1e5,
    const bool energy = true) {
  const double temperature = 0.25;
  const double beta_mu = -7.94;
  MonteCarlo mc;
  mc.set(MakeRandomMT19937({{"seed", "time"}}));
  mc.set(rpm({
    {"delta", "0.3"},
    {"charge_ratio", "2"},
    {"cubic_box_length", "15"},
    {"cutoff", "7.5"},
    {"kmax_squared", "25"},
    {"alpha", str(5./15)}}));
  if (min > 0) {
    ASSERT(min == 1, "unrecognized min: " << min);
    mc.get_system()->get_configuration()->add_particle_of_type(1);
  }
  mc.set(MakeThermoParams({{"beta", str(1/temperature)},
     {"chemical_potential0", str(beta_mu*temperature)},
     {"chemical_potential1", str(beta_mu*temperature)}}));
  auto criteria = MakeFlatHistogram(
    MakeMacrostateNumParticles(
      Histogram({{"width", "1"}, {"max", str(max)}, {"min", str(min)}})),
    MakeTransitionMatrix({{"min_sweeps", "100"}}),
    MakeAHalfB({{"extra", "1"}}));
  mc.set(criteria);
  mc.add(MakeCriteriaUpdater({{"steps_per", str(steps_per)}}));
  mc.add(MakeCriteriaWriter({
    {"steps_per", str(steps_per)},
    {"file_name", "tmp/dival_egce_crit.txt"},
  }));
  mc.add(MakeLogAndMovie({{"steps_per", str(steps_per)}, {"file_name", "tmp/dival_egce"}}));
  mc.add(MakeCheckProperties({{"steps_per", str(steps_per)}, {"tolerance", str(1e-12)}}));
  mc.add(MakeCheckEnergyAndTune({{"steps_per", str(steps_per)}, {"tolerance", str(1e-4)}}));
  const double charge_minus = mc.configuration().model_params().charge().value(1);
  mc.add(MakeCheckNetCharge({{"steps_per", str(steps_per)},
                             {"maximum", str(-charge_minus)},
                             {"minimum", str(charge_minus)}}));
  if (energy) {
    mc.add(MakeEnergy({
      {"file_name", "tmp/dival_egce_energy"},
      {"steps_per_update", "1"},
      {"steps_per_write", str(steps_per)},
      {"multistate", "true"}}));
  }
  return mc;
}

void compare_lnpi(const MonteCarlo& mc, const int min) {
  int shift = 0;
  if (min == 1) shift = -1;
  const LnProbability lnpi =
    FlatHistogram(mc.criteria()).bias().ln_prob().reduce(3, shift);
  INFO(feasst_str(lnpi.values()));
  int index = 0;
  if (min != 1) {
    EXPECT_NEAR(lnpi.value(index), -6.6615, 0.1);
    ++index;
  }
  EXPECT_NEAR(lnpi.value(index), -3.6256, 0.1);
  ++index;
  EXPECT_NEAR(lnpi.value(index), -2.1046, 0.1);
  ++index;
  EXPECT_NEAR(lnpi.value(index), -1.3685, 0.1);
  ++index;
  EXPECT_NEAR(lnpi.value(index), -1.1371, 0.1);
  ++index;
  EXPECT_NEAR(lnpi.value(index), -1.2911, 0.1);
}

void compare_lnpi_en(const MonteCarlo& mc, const int min) {
  compare_lnpi(mc, min);
  const std::vector<std::shared_ptr<Analyze> >& en =
    mc.analyzers().back()->analyzers();
  int index = 0;
  if (min != 1) {
    EXPECT_NEAR(en[0]->accumulator().average(), 0, 1e-14);
    ++index;
  }
  EXPECT_NEAR(en[index]->accumulator().average(), -0.0903901, 1e-4);
  ++index;
  EXPECT_NEAR(en[index]->accumulator().average(), -0.958108, 0.03);
  ++index;
  EXPECT_NEAR(en[index]->accumulator().average(), -1.30701, 0.03);
  ++index;
  EXPECT_NEAR(en[index]->accumulator().average(), -1.57156, 0.04);
  ++index;
  EXPECT_NEAR(en[index]->accumulator().average(), -2.60241, 0.05);
  ++index;
  EXPECT_NEAR(en[index]->accumulator().average(), -2.98115, 0.05);
  ++index;
  EXPECT_NEAR(en[index]->accumulator().average(), -3.30761, 0.06);
  ++index;
  EXPECT_NEAR(en[index]->accumulator().average(), -4.45619, 0.07);
  ++index;
  EXPECT_NEAR(en[index]->accumulator().average(), -4.85254, 0.08);
  ++index;
  EXPECT_NEAR(en[index]->accumulator().average(), -5.19956, 0.09);
  ++index;
  EXPECT_NEAR(en[index]->accumulator().average(), -6.38971, 0.12);
  ++index;
  EXPECT_NEAR(en[index]->accumulator().average(), -6.80956, 0.12);
  ++index;
  EXPECT_NEAR(en[index]->accumulator().average(), -7.18859, 0.12);
  ++index;
  EXPECT_NEAR(en[index]->accumulator().average(), -8.43124, 0.15);
  ++index;
  EXPECT_NEAR(en[index]->accumulator().average(), -8.85025, 0.16);
}

TEST(MonteCarlo, rpm_egce_divalent_LONG) {
  const int min = 0;
  MonteCarlo mc = dival_egce(min);
  mc.add(MakeTrialTranslate({{"weight", "0.25"}, {"tunable_param", "0.1"}}));
  mc.add(MakeTrialTransfer({{"weight", "1."}, {"particle_type", "0"}}));
  mc.add(MakeTrialTransfer({{"weight", "1."}, {"particle_type", "1"}}));
  mc.run_until_complete();
  compare_lnpi_en(mc, min);
}

TEST(MonteCarlo, rpm_egce_min1_divalent_LONG) {
  const int min = 1;
  MonteCarlo mc = dival_egce(1);
  mc.add(MakeTrialTranslate({{"weight", "0.25"}, {"tunable_param", "0.1"}}));
  mc.add(MakeTrialTransfer({{"weight", "1."}, {"particle_type", "0"}}));
  mc.add(MakeTrialTransfer({{"weight", "1."}, {"particle_type", "1"}}));
  mc.run_until_complete();
  compare_lnpi_en(mc, min);
}

TEST(MonteCarlo, rpm_egce_avb_divalent_LONG) {
  const int min = 1;
  MonteCarlo mc = dival_egce(min);
  mc.add(MakeTrialTranslate({{"weight", "0.25"}, {"tunable_param", "0.1"}}));
  mc.set(MakeRandomMT19937({{"seed", "time"}}));
  mc.set(1, MakePotential(MakeModelTwoBodyFactory({MakeHardSphere(),
                                                   MakeChargeScreened()}),
                          MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
//  mc.set(MakeThermoParams({{"beta", str(mc.criteria().beta())},
//     {"chemical_potential0", str(mc.criteria().chemical_potential(0))},
//     {"chemical_potential1", str(mc.criteria().chemical_potential(1))}}));
  mc.set(MakeFlatHistogram(
    MakeMacrostateNumParticles(
      Histogram({{"width", "1"}, {"max", "15"}, {"min", str(min)}})),
    MakeTransitionMatrix({{"min_sweeps", "100"}}),
    MakeAHalfB({{"extra", "1"}})));
  mc.add(MakeNeighborCriteria({{"maximum_distance", "7.5"},
                               {"minimum_distance", "1"},
                               {"site_type0", "0"},
                               {"site_type1", "1"},
                               {"potential_index", "1"}}));
  mc.add(MakeTrialTransferAVB({
    {"neighbor_index", "0"},
    {"weight", "1."},
    {"particle_type", "0"},
    {"target_particle_type", "1"}}));
  mc.add(MakeTrialTransferAVB({
    {"neighbor_index", "0"},
    {"weight", "1."},
    {"particle_type", "1"},
    {"target_particle_type", "0"}}));
  mc.run_until_complete();
  compare_lnpi_en(mc, min);
}

TEST(MonteCarlo, rpm_egce_divalent_avb_and_not) {
  const int min = 0, max = 15, steps_per = 1e3;
  MonteCarlo mc = dival_egce(min, max, steps_per);
  mc.add(MakeTrialTranslate({{"weight", "0.25"}, {"tunable_param", "0.1"}}));
  mc.set(MakeRandomMT19937({{"seed", "time"}}));
  mc.set(1, MakePotential(MakeModelTwoBodyFactory({MakeHardSphere(),
                                                   MakeChargeScreened()}),
                          MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
  mc.add(MakeNeighborCriteria({{"maximum_distance", "7.5"},
                               {"minimum_distance", "1"},
                               {"site_type0", "0"},
                               {"site_type1", "1"},
                               {"potential_index", "1"}}));
  mc.add(MakeTrialTransfer({{"particle_type", "0"}}));
  mc.add(MakeTrialTransfer({{"particle_type", "1"}}));
  mc.add(MakeTrialTransferAVB({
    {"neighbor_index", "0"},
    {"weight", "1."},
    {"particle_type", "0"},
    {"target_particle_type", "1"}}));
  mc.add(MakeTrialTransferAVB({
    {"neighbor_index", "0"},
    {"weight", "1."},
    {"particle_type", "1"},
    {"target_particle_type", "0"}}));
  EXPECT_TRUE(mc.trial(7).stage(0).select().is_ghost());
  mc.attempt(2*steps_per);
}

TEST(MonteCarlo, rpm_divalent_avb_VERY_LONG) {
  const int min = 0, max = 5, steps_per = 1e3;
  MonteCarlo mc = dival_egce(min, max, steps_per);
  mc.add(MakeTrialTranslate({{"weight", "0.25"}, {"tunable_param", "0.1"}}));
  mc.set(MakeRandomMT19937({{"seed", "time"}}));
  mc.set(1, MakePotential(MakeModelTwoBodyFactory({MakeHardSphere(),
                                                   MakeChargeScreened()}),
                          MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
  mc.add_to_reference(MakePotential(MakeDontVisitModel()));
//  mc.set(MakeThermoParams({{"beta", str(mc.criteria().beta())},
//     {"chemical_potential0", str(mc.criteria().chemical_potential(0))},
//     {"chemical_potential1", str(mc.criteria().chemical_potential(1))}}));
  auto criteria = MakeFlatHistogram(
    MakeMacrostateNumParticles(
      Histogram({{"width", "1"}, {"max", str(max)}, {"min", str(min)}}),
      {{"particle_type", "0"}}),
    MakeTransitionMatrix({{"min_sweeps", "1000"}}));
  mc.set(criteria);
  mc.add(MakeTrialTranslate({{"weight", "0.25"}, {"tunable_param", "0.1"}}));
  mc.add(MakeNeighborCriteria({{"maximum_distance", "4"},
                               {"minimum_distance", "1"},
                               {"site_type0", "0"},
                               {"site_type1", "1"},
                               {"potential_index", "1"}}));
  mc.add(MakeTrialTransferAVBDivalent({
    {"neighbor_index", "0"},
    {"weight", "1."},
    {"particle_type", "0"},
    {"particle_type_a", "1"},
    {"particle_type_b", "1"},
    {"reference_index", "0"}}));
  mc.run_until_complete();
  const LnProbability lnpi = FlatHistogram(mc.criteria()).bias().ln_prob();
  EXPECT_NEAR(lnpi.value(0), -6.7005955776549158, 0.09);
  EXPECT_NEAR(lnpi.value(1), -3.6523345299136007, 0.06);
  EXPECT_NEAR(lnpi.value(2), -2.1178631459398805, 0.03);
  EXPECT_NEAR(lnpi.value(3), -1.3652342629553453, 0.02);
  EXPECT_NEAR(lnpi.value(4), -1.1336431696116527, 0.02);
  EXPECT_NEAR(lnpi.value(5), -1.289634124762612, 0.02);
  const std::vector<std::shared_ptr<Analyze> >& en =
    mc.analyzers().back()->analyzers();
  EXPECT_NEAR(en[0]->accumulator().average(), 0, 1e-14);
  EXPECT_NEAR(en[1]->accumulator().average(), -1.3278876302141585, 0.05);
  EXPECT_NEAR(en[2]->accumulator().average(), -3.0162868737745732, 0.05);
  EXPECT_NEAR(en[3]->accumulator().average(), -4.8648645814174927, 0.05);
  EXPECT_NEAR(en[4]->accumulator().average(), -6.8089768188067694, 0.05);
  EXPECT_NEAR(en[5]->accumulator().average(), -8.8377616317395002, 0.05);
}

TEST(MonteCarlo, rpm_divalent_morph_LONG) {
  const int min = 0, max = 5, steps_per = 1e3;
  MonteCarlo mc = dival_egce(min, max, steps_per, false);
  mc.set(MakeRandomMT19937({{"seed", "123"}}));

  // add growth expanded particle types with half charge
  { Configuration * config = mc.get_system()->get_configuration();
    const double q_plus = config->model_params().charge().value(0);
    const double q_minus = config->model_params().charge().value(1);
    config->add_particle_type(install_dir() +
      "/plugin/charge/forcefield/rpm_plus.fstprt", "0.5");
    config->add_particle_type(install_dir() +
      "/plugin/charge/forcefield/rpm_minus.fstprt", "0.5");
    config->set_model_param("charge", 2, 0.5*q_plus);
    config->set_model_param("charge", 3, 0.5*q_minus);
    config->set_model_param("cutoff", 2, 7.5);
    config->set_model_param("cutoff", 3, 7.5);
    //mc.get_system()->precompute();  // update charge
  }
  mc.add_to_reference(MakePotential(MakeDontVisitModel()));
  const std::vector<std::vector<int> > grow_sequence = {{2, 3, 3}, {0, 1, 1}};
  mc.set(MakeThermoParams({{"beta", str(mc.system().thermo_params().beta())},
    {"chemical_potential0", str(mc.system().thermo_params().chemical_potential(0))},
    {"chemical_potential1", str(mc.system().thermo_params().chemical_potential(1))},
    {"chemical_potential2", str(mc.system().thermo_params().chemical_potential(0))},
    {"chemical_potential3", str(mc.system().thermo_params().chemical_potential(1))}}));
  mc.set(MakeFlatHistogram(
    MakeMacrostateMorph(
      grow_sequence,
      Histogram({{"width", str(1./grow_sequence.size())},
                 {"max", str(max)}, {"min", str(min)}})),
    // MakeWangLandau({{"min_flatness", "25"}}),
    MakeTransitionMatrix({{"min_sweeps", "1000"}})));
  mc.initialize_criteria();
  mc.initialize_analyzers();
  mc.add(MakeTrialTranslate({{"weight", "0.25"}, {"tunable_param", "0.1"}}));
  mc.add(MakeTrialMorphExpanded(grow_sequence, {{"reference_index", "0"}}));
  mc.add(MakeEnergy({
    {"file_name", "tmp/dival_egce_energy"},
    {"steps_per_update", "1"},
    {"steps_per_write", str(steps_per)},
    {"multistate", "true"}}));
  mc.run_until_complete();
  const LnProbability lnpi = FlatHistogram(mc.criteria()).bias().ln_prob().reduce(2);
//  EXPECT_NEAR(lnpi.value(0), -6.6615, 0.1);
//  EXPECT_NEAR(lnpi.value(1), -3.6256, 0.1);
//  EXPECT_NEAR(lnpi.value(2), -2.02966, 0.1);
//  EXPECT_NEAR(lnpi.value(3), -1.35573, 0.1);
//  EXPECT_NEAR(lnpi.value(4), -1.16195, 0.1);
//  EXPECT_NEAR(lnpi.value(5), -1.34209, 0.1);
  EXPECT_NEAR(lnpi.value(0), -6.7005955776549158, 0.09);
  EXPECT_NEAR(lnpi.value(1), -3.6523345299136007, 0.06);
  EXPECT_NEAR(lnpi.value(2), -2.1178631459398805, 0.03);
  EXPECT_NEAR(lnpi.value(3), -1.3652342629553453, 0.03);
  EXPECT_NEAR(lnpi.value(4), -1.1336431696116527, 0.025);
  EXPECT_NEAR(lnpi.value(5), -1.289634124762612, 0.025);
  const std::vector<std::shared_ptr<Analyze> >& en =
    mc.analyzers().back()->analyzers();
//  EXPECT_NEAR(en[0]->accumulator().average(), 0, 1e-14);
//  EXPECT_NEAR(en[2]->accumulator().average(), -1.30701, 0.05);
//  EXPECT_NEAR(en[4]->accumulator().average(), -2.98115, 0.05);
//  EXPECT_NEAR(en[6]->accumulator().average(), -4.85254, 0.08);
//  EXPECT_NEAR(en[8]->accumulator().average(), -6.80956, 0.12);
//  EXPECT_NEAR(en[10]->accumulator().average(), -8.85025, 0.16);
  EXPECT_NEAR(en[0]->accumulator().average(), 0, 1e-14);
  EXPECT_NEAR(en[2]->accumulator().average(), -1.3278876302141585, 0.05);
  EXPECT_NEAR(en[4]->accumulator().average(), -3.0162868737745732, 0.05);
  EXPECT_NEAR(en[6]->accumulator().average(), -4.8648645814174927, 0.05);
  EXPECT_NEAR(en[8]->accumulator().average(), -6.8089768188067694, 0.05);
  EXPECT_NEAR(en[10]->accumulator().average(), -8.8377616317395002, 0.06);
}

double energy_av467(const int macro, const MonteCarlo& mc) {
  return mc.analyzers().back()->analyzers()[macro]->accumulator().average();
}

//TEST(MonteCarlo, lj_fh_trial_grow_LONG) {
//  MonteCarlo mc;
//  mc.set(MakeRandomMT19937({{"seed", "123"}}));
//  mc.set(lennnnnard_jones({{"dual_cut", "1"}}));
//  mc.get_system()->get_configuration()->add_particle_of_type(0);
//  mc.set(MakeThermoParams({{"beta", str(1/1.5)}, {"chemical_potential", "-2.352321"}}));
//  mc.set(MakeFlatHistogram(
//    MakeMacrostateNumParticles(Histogram({{"width", "1"}, {"max", "5"}, {"min", "1"}})),
//    MakeTransitionMatrix({{"min_sweeps", "1000"}})));
//  const std::string ns = "4";
//  const std::string ref = "0";
//  //const std::string ns = "1";
//  //const std::string ref = "-1";
//  mc.add(MakeTrialGrow({{{"particle_type", "0"}, {"translate", "true"}, {"tunable_param", "1"}, {"site", "0"}, {"num_steps", ns}, {"reference_index", ref}}}));
//  mc.add(MakeTrialGrow({{{"particle_type", "0"}, {"transfer", "true"}, {"weight", "4"}, {"site", "0"}, {"num_steps", ns}, {"reference_index", ref}}}));
//  //mc.add(MakeTrialTranslate());
//  //mc.add(MakeTrialTransfer({{"particle_type", "0"}, {"weight", "4"}}));
//  mc.add(MakeCheckEnergyAndTune({{"steps_per", "1e0"}}));
//  mc.add(MakeCriteriaUpdater({{"steps_per", "1e3"}}));
//  mc.add(MakeCriteriaWriter({{"steps_per", "1e3"}, {"file_name", "tmp/ljcrit.txt"}}));
//  mc.add(MakeEnergy({{"file_name", "tmp/lj_fh_energy"},
//    {"steps_per_update", "1"},
//    {"steps_per_write", "1e3"},
//    {"multistate", "true"}}));
//  mc.run_until_complete();
//  const LnProbability lnpi = FlatHistogram(mc.criteria()).bias().ln_prob();
//  EXPECT_NEAR(lnpi.value(0), -14.037373358321800000, 0.02);
//  EXPECT_NEAR(lnpi.value(1), -10.050312091655200000, 0.02);
//  EXPECT_NEAR(lnpi.value(2), -6.458920624988570000, 0.02);
//  EXPECT_NEAR(lnpi.value(3), -3.145637424988510000, 0.01);
//  EXPECT_NEAR(lnpi.value(4), -0.045677458321876000, 0.005);
//  EXPECT_NEAR(energy_av467(0, mc), -0.000605740233333333, 1e-8);
//  EXPECT_NEAR(energy_av467(1, mc), -0.030574223333333334, 0.001);
//  EXPECT_NEAR(energy_av467(2, mc), -0.089928316, 0.002);
//  EXPECT_NEAR(energy_av467(3, mc), -0.1784570533333333, 0.004);
//  EXPECT_NEAR(energy_av467(4, mc), -0.29619201333333334, 0.006);
//}

TEST(MonteCarlo, lj_fh_trial_grow_liquid_LONG) {
  MonteCarlo mc;
  mc.add(MakeConfiguration({{"cubic_box_length", "8"}, {"particle_type0", "../forcefield/lj.fstprt"}}));
  mc.add(MakePotential(MakeLennardJones()));
  mc.add(MakePotential(MakeLongRangeCorrections()));
  mc.run(MakeAddReference({{"cutoff", "1"}, {"use_cell", "true"}}));
  mc.get_system()->get_configuration()->add_particle_of_type(0);
  mc.set(MakeThermoParams({{"beta", str(1/1.5)}, {"chemical_potential", "-2.352321"}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialAdd({{"particle_type", "0"}}));
  mc.run(MakeRun({{"until_num_particles", "100"}}));
  mc.run(MakeRemoveTrial({{"name", "TrialAdd"}}));
  mc.set(MakeFlatHistogram(
    MakeMacrostateNumParticles(Histogram({{"width", "1"}, {"max", "105"}, {"min", "100"}})),
    MakeTransitionMatrix({{"min_sweeps", "1000"}})));
  const std::string ns = "4";
  const std::string ref = "0";
//  const std::string ns = "1";
//  const std::string ref = "-1";
  mc.add(MakeTrialGrow({{{"particle_type", "0"}, {"translate", "true"}, {"tunable_param", "1"}, {"site", "0"}, {"num_steps", ns}, {"reference_index", ref}}}));
  mc.add(MakeTrialTransfer({{"particle_type", "0"}, {"weight", "4"}}));
  //mc.add(MakeTrialGrow({{{"particle_type", "0"}, {"transfer", "true"}, {"weight", "4"}, {"site", "0"}, {"num_steps", ns}, {"reference_index", ref}}}));
  //mc.add(MakeTrialTranslate());
  const std::string steps_per = "1e3";
  mc.add(MakeCheckEnergyAndTune({{"steps_per", steps_per}}));
  mc.add(MakeLogAndMovie({{"steps_per", steps_per}, {"file_name", "tmp/lj"}}));
  mc.add(MakeCriteriaUpdater({{"steps_per", steps_per}}));
  mc.add(MakeCriteriaWriter({{"steps_per", steps_per}, {"file_name", "tmp/ljcrit.txt"}}));
  mc.add(MakeEnergy({{"file_name", "tmp/lj_fh_energy"},
    {"steps_per_update", "1"},
    {"steps_per_write", steps_per},
    {"multistate", "true"}}));
  mc.run_until_complete();
  const LnProbability lnpi = FlatHistogram(mc.criteria()).bias().ln_prob();
  EXPECT_NEAR(lnpi.value(0), -4.92194963175925, 0.025);
  EXPECT_NEAR(lnpi.value(1), -4.03855513175926, 0.02);
  EXPECT_NEAR(lnpi.value(2), -3.15822813175925, 0.02);
  EXPECT_NEAR(lnpi.value(3), -2.28019483175925, 0.015);
  EXPECT_NEAR(lnpi.value(4), -1.40647303175926, 0.005);
  EXPECT_NEAR(lnpi.value(5), -0.535594831759248, 0.005);
  EXPECT_NEAR(energy_av467(0, mc), -1.381223800E+02, 0.5);
  EXPECT_NEAR(energy_av467(1, mc), -1.408257000E+02, 0.5);
  EXPECT_NEAR(energy_av467(2, mc), -1.435426000E+02, 0.5);
  EXPECT_NEAR(energy_av467(3, mc), -1.462802100E+02, 0.5);
  EXPECT_NEAR(energy_av467(4, mc), -1.490501400E+02, 0.5);
  EXPECT_NEAR(energy_av467(5, mc), -1.518517300E+02, 0.6);
}

}  // namespace feasst
