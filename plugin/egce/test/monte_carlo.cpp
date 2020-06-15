#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "system/include/hard_sphere.h"
#include "system/include/model_two_body_factory.h"
#include "system/include/dont_visit_model.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/trial_rotate.h"
#include "monte_carlo/include/trial_transfer.h"
#include "monte_carlo/include/trial_add.h"
#include "monte_carlo/include/seek_num_particles.h"
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
#include "ewald/include/check_net_charge.h"
#include "ewald/include/utils.h"
#include "ewald/include/charge_screened.h"
#include "cluster/include/energy_map_all.h"
#include "cluster/include/trial_transfer_avb.h"
#include "cluster/include/trial_transfer_avb_divalent.h"
#include "egce/include/a_equal_b.h"
#include "egce/include/a_half_b.h"

namespace feasst {

MonteCarlo rpm_egce(const int min = 0, const std::string dual_cut = "-1") {
  MonteCarlo mc;
  // mc.set(MakeRandomMT19937({{"seed", "default"}}));
  mc.set(rpm({
    {"cubic_box_length", "12"},
    {"cutoff", "4.891304347826090"},
    {"alphaL", "6.87098396396261"},
    {"dual_cut", dual_cut}}));
  const double temperature = 0.047899460618081;
  const double beta_mu = -13.94;
  if (min > 0) {
    mc.set(MakeMetropolis({{"beta", "0.01"}, {"chemical_potential", "1"}}));
    SeekNumParticles(min).with_trial_add().run(&mc);
  }
  auto criteria = MakeFlatHistogram(
    MakeMacrostateNumParticles(
      Histogram({{"width", "1"}, {"max", "4"}, {"min", str(min)}})),
    MakeTransitionMatrix({{"min_sweeps", "100"}}),
    MakeAEqualB({{"extra_A", "1"}}),
    { {"beta", str(1/temperature)},
      {"chemical_potential0", str(beta_mu*temperature)},
      {"chemical_potential1", str(beta_mu*temperature)}});
  mc.set(criteria);
  const int steps_per = 1e5;
  mc.add(MakeCriteriaUpdater({{"steps_per", str(steps_per)}}));
  mc.add(MakeCriteriaWriter({{"steps_per", str(steps_per)},
                             {"file_name", "tmp/rpm_egce_crit.txt"}}));
  mc.add(MakeLogAndMovie({{"steps_per", str(steps_per)}, {"file_name", "tmp/rpm_egce"}}));
  mc.add(MakeCheckEnergyAndTune({{"steps_per", str(steps_per)}, {"tolerance", "1e-8"}}));
  // mc.add(MakeCheckProperties({{"steps_per", str(steps_per)}}));
  // mc.add(MakeCPUTime({{"steps_per", str(5*steps_per)}}));
  mc.add(MakeCheckNetCharge({{"maximum", "1."}, {"minimum", str(-NEAR_ZERO)}}));
  mc.add(MakeEnergy({
    {"file_name", "tmp/rpm_egce_energy"},
    {"steps_per_update", "1"},
    {"steps_per_write", str(steps_per)},
    {"multistate", "true"}}));
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
    mc.add(MakeTrialTranslate({
      {"weight", "0.25"},
      {"tunable_param", "0.1"},
      {"reference_index", ref},
      {"num_steps", str(num_steps)}}));
    mc.add(MakeTrialTransfer({
      {"weight", "1."},
      {"particle_type", "0"},
      {"reference_index", ref},
      {"num_steps", str(num_steps)}}));
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
    EXPECT_NEAR(lnpi3.value(0), -1.2994315780357, 0.1);
    EXPECT_NEAR(lnpi3.value(1), -1.08646312498868, 0.1);
    EXPECT_NEAR(lnpi3.value(2), -0.941850889679828, 0.06);
    const std::vector<std::shared_ptr<Analyze> >& en = mc2.analyzers().back()->analyzers();
    EXPECT_NEAR(en[0]->accumulator().average(), 0, 1e-14);
    EXPECT_NEAR(en[1]->accumulator().average(), -0.115474, 1e-6);
    EXPECT_NEAR(en[2]->accumulator().average(), -0.939408, 0.02);
    EXPECT_NEAR(en[3]->accumulator().average(), -1.32485, 0.03);
    EXPECT_NEAR(en[4]->accumulator().average(), -2.02625, 0.04);
  }
}

TEST(MonteCarlo, rpm_egce_fh_min1_LONG) {
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
  EXPECT_NEAR(lnpi.value(0), -5.05, 0.1);
  EXPECT_NEAR(lnpi.value(1), -0.77327, 0.1);
  EXPECT_NEAR(lnpi.value(2), -3.55107, 0.06);
  EXPECT_NEAR(lnpi.value(3), -0.686417, 0.06);
  const std::vector<std::shared_ptr<Analyze> >& en = mc.analyzers().back()->analyzers();
  EXPECT_NEAR(en[0]->accumulator().average(), -0.115474, 1e-6);
  EXPECT_NEAR(en[1]->accumulator().average(), -0.939408, 0.02);
  EXPECT_NEAR(en[2]->accumulator().average(), -1.32485, 0.03);
  EXPECT_NEAR(en[3]->accumulator().average(), -2.02625, 0.04);
}

TEST(MonteCarlo, rpm_egce_avb_fh_LONG) {
  MonteCarlo mc = rpm_egce(1);
  mc.add(MakeTrialTranslate({{"weight", "0.25"}, {"tunable_param", "0.1"}}));
  // mc.set(MakeRandomMT19937({{"seed", "1346867550"}}));
  mc.set(1, Potential(MakeModelTwoBodyFactory({MakeHardSphere(),
                                               MakeChargeScreened()}),
                      MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
  mc.get_system()->energy();
  auto neighbor_criteria = MakeNeighborCriteria({{"maximum_distance", "3"},
                                                 {"minimum_distance", "1"},
                                                 {"site_type0", "0"},
                                                 {"site_type1", "1"},
                                                 {"potential_index", "1"}});
  mc.add(MakeTrialTransferAVB(neighbor_criteria,
    { {"weight", "1."},
      {"particle_type", "0"},
      {"target_particle_type", "1"}}));
  mc.add(MakeTrialTransferAVB(neighbor_criteria,
    { {"weight", "1."},
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
  EXPECT_NEAR(en[3]->accumulator().average(), -2.02625, 0.04);
}

MonteCarlo dival_egce(
    const int min = 0,
    const int max = 15,
    const int steps_per = 1e5) {
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
    {"alphaL", "5"}}));
  if (min > 0) {
    ASSERT(min == 1, "unrecognized min: " << min);
    mc.get_system()->get_configuration()->add_particle_of_type(1);
  }
  auto criteria = MakeFlatHistogram(
    MakeMacrostateNumParticles(
      Histogram({{"width", "1"}, {"max", str(max)}, {"min", str(min)}})),
    MakeTransitionMatrix({{"min_sweeps", "100"}}),
    MakeAHalfB({{"extra", "1"}}),
    {{"beta", str(1/temperature)},
     {"chemical_potential0", str(beta_mu*temperature)},
     {"chemical_potential1", str(beta_mu*temperature)}});
  mc.set(criteria);
  mc.add(MakeCriteriaUpdater({{"steps_per", str(steps_per)}}));
  mc.add(MakeCriteriaWriter({
    {"steps_per", str(steps_per)},
    {"file_name", "tmp/dival_egce_crit.txt"},
  }));
  mc.add(MakeLogAndMovie({{"steps_per", str(steps_per)}, {"file_name", "tmp/dival_egce"}}));
  mc.add(MakeCheckEnergyAndTune({{"steps_per", str(steps_per)}, {"tolerance", str(1e-4)}}));
  const double charge_minus = mc.configuration().model_params().charge().value(1);
  mc.add(MakeCheckNetCharge({{"steps_per", str(steps_per)},
                             {"maximum", str(-charge_minus)},
                             {"minimum", str(charge_minus)}}));
  mc.add(MakeEnergy({
    {"file_name", "tmp/dival_egce_energy"},
    {"steps_per_update", "1"},
    {"steps_per_write", str(steps_per)},
    {"multistate", "true"}}));
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
  EXPECT_NEAR(en[index]->accumulator().average(), -1.57156, 0.03);
  ++index;
  EXPECT_NEAR(en[index]->accumulator().average(), -2.60241, 0.04);
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
  mc.set(1, Potential(MakeModelTwoBodyFactory({MakeHardSphere(),
                                               MakeChargeScreened()}),
                      MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
  mc.set(MakeFlatHistogram(
    MakeMacrostateNumParticles(
      Histogram({{"width", "1"}, {"max", "15"}, {"min", str(min)}})),
    MakeTransitionMatrix({{"min_sweeps", "100"}}),
    MakeAHalfB({{"extra", "1"}}),
    {{"beta", str(mc.criteria().beta())},
     {"chemical_potential0", str(mc.criteria().chemical_potential(0))},
     {"chemical_potential1", str(mc.criteria().chemical_potential(1))}}));
  auto neighbor_criteria = MakeNeighborCriteria({{"maximum_distance", "7.5"},
                                                 {"minimum_distance", "1"},
                                                 {"site_type0", "0"},
                                                 {"site_type1", "1"},
                                                 {"potential_index", "1"}});
  mc.add(MakeTrialTransferAVB(neighbor_criteria,
    { {"weight", "1."},
      {"particle_type", "0"},
      {"target_particle_type", "1"}}));
  mc.add(MakeTrialTransferAVB(neighbor_criteria,
    { {"weight", "1."},
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
  mc.set(1, Potential(MakeModelTwoBodyFactory({MakeHardSphere(),
                                               MakeChargeScreened()}),
                      MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
  auto neighbor_criteria = MakeNeighborCriteria({{"maximum_distance", "7.5"},
                                                 {"minimum_distance", "1"},
                                                 {"site_type0", "0"},
                                                 {"site_type1", "1"},
                                                 {"potential_index", "1"}});
  mc.add(MakeTrialTransfer({{"particle_type", "0"}}));
  mc.add(MakeTrialTransfer({{"particle_type", "1"}}));
  mc.add(MakeTrialTransferAVB(neighbor_criteria,
    { {"weight", "1."},
      {"particle_type", "0"},
      {"target_particle_type", "1"}}));
  mc.add(MakeTrialTransferAVB(neighbor_criteria,
    { {"weight", "1."},
      {"particle_type", "1"},
      {"target_particle_type", "0"}}));
  mc.attempt(2*steps_per);
}

TEST(MonteCarlo, rpm_divalent_avb_LONG) {
  const int min = 0, max = 5, steps_per = 1e3;
  MonteCarlo mc = dival_egce(min, max, steps_per);
  mc.add(MakeTrialTranslate({{"weight", "0.25"}, {"tunable_param", "0.1"}}));
  mc.set(MakeRandomMT19937({{"seed", "time"}}));
  mc.set(1, Potential(MakeModelTwoBodyFactory({MakeHardSphere(),
                                               MakeChargeScreened()}),
                      MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
  mc.add_to_reference(Potential(MakeDontVisitModel()));
    SeekNumParticles(3*max)
      .with_metropolis(
        MakeAHalfB({{"extra", "1"}}),
        { {"beta", "0.1"},
          {"chemical_potential0", "1"},
          {"chemical_potential1", "1"}})
      .add(MakeTrialAdd({{"particle_type", "0"}}))
      .add(MakeTrialAdd({{"particle_type", "1"}}))
      .run(&mc);
  auto criteria = MakeFlatHistogram(
    MakeMacrostateNumParticles(
      Histogram({{"width", "1"}, {"max", str(max)}, {"min", str(min)}}),
      {{"particle_type", "0"}}),
    MakeTransitionMatrix({{"min_sweeps", "100"}}),
    {{"beta", str(mc.criteria().beta())},
     {"chemical_potential0", str(mc.criteria().chemical_potential(0))},
     {"chemical_potential1", str(mc.criteria().chemical_potential(1))}});
  mc.set(criteria);
  auto neighbor_criteria = MakeNeighborCriteria({{"maximum_distance", "1.5"},
                                                 {"minimum_distance", "1"},
                                                 {"site_type0", "0"},
                                                 {"site_type1", "1"},
                                                 {"potential_index", "1"}});
  mc.add(MakeTrialTransferAVBDivalent(neighbor_criteria,
    { {"weight", "1."},
      {"particle_type", "0"},
      {"particle_type_a", "1"},
      {"particle_type_b", "1"},
      {"reference_index", "0"}}));
  mc.run_until_complete();
  const LnProbability lnpi = FlatHistogram(mc.criteria()).bias().ln_prob();
  EXPECT_NEAR(lnpi.value(0), -5.41, 0.1);
  EXPECT_NEAR(lnpi.value(1), -3.42548, 0.1);
  EXPECT_NEAR(lnpi.value(2), -2.02966, 0.1);
  EXPECT_NEAR(lnpi.value(3), -1.35573, 0.1);
  EXPECT_NEAR(lnpi.value(4), -1.16195, 0.1);
  EXPECT_NEAR(lnpi.value(5), -1.34209, 0.1);
  const std::vector<std::shared_ptr<Analyze> >& en =
    mc.analyzers().back()->analyzers();
  EXPECT_NEAR(en[0]->accumulator().average(), 0, 1e-14);
  EXPECT_NEAR(en[1]->accumulator().average(), -1.462473, 0.03);
  EXPECT_NEAR(en[2]->accumulator().average(), -3.015735, 0.05);
  EXPECT_NEAR(en[3]->accumulator().average(), -4.879190, 0.08);
  EXPECT_NEAR(en[4]->accumulator().average(), -6.829874, 0.12);
  EXPECT_NEAR(en[5]->accumulator().average(), -8.815852, 0.16);
}

}  // namespace feasst
