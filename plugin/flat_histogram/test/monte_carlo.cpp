#include "utils/test/utils.h"
#include "math/include/histogram.h"
#include "math/include/random_mt19937.h"
#include "system/include/hard_sphere.h"
#include "system/include/dont_visit_model.h"
#include "system/include/model_two_body_factory.h"
#include "system/include/utils.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/test/monte_carlo_test.h"
#include "monte_carlo/include/seek_num_particles.h"
#include "monte_carlo/include/trial_rotate.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/trial_transfer.h"
#include "steppers/include/energy.h"
#include "steppers/include/criteria_writer.h"
#include "steppers/include/criteria_updater.h"
#include "steppers/include/cpu_time.h"
#include "steppers/include/check_properties.h"
#include "steppers/include/check_physicality.h"
#include "steppers/include/check_energy_and_tune.h"
#include "steppers/include/log_and_movie.h"
#include "ewald/include/ewald.h"
#include "ewald/include/charge_screened.h"
#include "ewald/include/charge_self.h"
#include "ewald/include/check_net_charge.h"
#include "ewald/include/utils.h"
#include "ewald/include/trial_transfer_multiple.h"
#include "flat_histogram/include/flat_histogram.h"
#include "flat_histogram/include/macrostate_num_particles.h"
#include "flat_histogram/include/transition_matrix.h"
#include "flat_histogram/include/wang_landau.h"
#include "flat_histogram/include/wltm.h"

namespace feasst {

double energy_av(const int macro, const MonteCarlo& mc) {
  return mc.analyzers().back()->analyzers()[macro]->accumulator().average();
}

//TEST(FlatHistogram, order) {
//  auto criteria = MakeFlatHistogram({
//    {"beta", str(1./1.5)},
//    {"chemical_potential", "-2.352321"}
//  });
////  TRY(
////    auto criteria2 = criteria;
////    criteria2->set(MakeWangLandau({{"min_flatness", "20"}}));
////    CATCH_PHRASE("set macrostate before bias");
////  );
//}

TEST(MonteCarlo, lj_fh) {
  for (int num_steps : {1, 2}) {
    for (const std::string bias_name : {"TM", "WL", "WLTM"}) {
      MonteCarlo mc;
      // mc.set(MakeRandomMT19937({{"seed", "default"}}));
      int ref = -1;
      if (num_steps == 1) {
        mc.set(lennard_jones());
      } else {
        ref = 0;
        mc.set(lennard_jones({{"dual_cut", "1."}}));
      }
      mc.set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "1."}}));
      mc.add(MakeTrialTranslate({
        {"weight", "1."},
        {"tunable_param", "1."},
        {"reference_index", str(ref)},
        {"num_steps", str(num_steps)}}));
      SeekNumParticles(1).with_trial_add().run(&mc);
      mc.add(MakeTrialTransfer({
        {"particle_type", "0"},
        {"reference_index", str(ref)},
        {"num_steps", str(num_steps)},
        {"weight", "4"}}));
      std::shared_ptr<Bias> bias;
      if (bias_name == "TM") {
        bias = MakeTransitionMatrix({{"min_sweeps", "10"}});
      } else if (bias_name == "WL") {
        bias = MakeWangLandau({{"min_flatness", "20"}});
      } else if (bias_name == "WLTM") {
        bias = MakeWLTM({{"collect_flatness", "15"},
          {"min_flatness", "20"}, {"min_sweeps", "10"}});
      } else {
        FATAL("unrecognized");
      }
      auto criteria = MakeFlatHistogram(
        MakeMacrostateNumParticles(
          Histogram({{"width", "1"}, {"max", "5"}, {"min", "1"}})),
        bias,
        { {"beta", str(1./1.5)},
          {"chemical_potential", "-2.352321"}});
  //      {{"soft_max", "5"}, {"soft_min", "1"}}));
  //      {{"particle_type", "0"}}));
      INFO(criteria->bias().class_name());
      mc.set(criteria);
      const std::string steps_per(str(1e4));
      // const std::string steps_per(str(1e4));
      mc.add(MakeLogAndMovie({{"steps_per", str(steps_per)}, {"file_name", "tmp/lj_fh"}}));
      mc.add(MakeCheckEnergyAndTune({{"steps_per", str(steps_per)}}));
      mc.add(MakeCriteriaUpdater({{"steps_per", str(1)}}));
      mc.add(MakeCriteriaWriter({
        {"steps_per", steps_per},
        {"file_name", "tmp/ljcrit.txt"}}));
      auto energy = MakeEnergy({
        {"file_name", "tmp/lj_fh_energy"},
        {"steps_per_update", "1"},
        {"steps_per_write", str(steps_per)},
        {"multistate", "true"}});
      EXPECT_EQ(energy->steps_per_update(), 1);
      EXPECT_EQ(energy->steps_per_write(), 1e4);
      mc.add(energy);
      //mc.attempt(1e4);
      mc.attempt(1e5); // note more than 1e4 steps required for TM
      //mc.run_until_complete();
      // INFO(mc.criteria().write());

      //MonteCarlo mc2 = test_serialize_no_comp(mc);
      test_serialize(mc);

      // compare with known values of lnpi
      const LnProbability * lnpi = &criteria->bias().ln_prob();
      //EXPECT_NEAR(lnpi->value(0), -18.707570324988800000, 0.55);
      EXPECT_NEAR(lnpi->value(0), -14.037373358321800000, 0.6);
      EXPECT_NEAR(lnpi->value(1), -10.050312091655200000, 0.6);
      EXPECT_NEAR(lnpi->value(2), -6.458920624988570000, 0.55);
      EXPECT_NEAR(lnpi->value(3), -3.145637424988510000, 0.55);
      EXPECT_NEAR(lnpi->value(4), -0.045677458321876000, 0.55);

      // compare with known values of energy
      //EXPECT_NEAR(energy_av(0, mc), 0, 1e-14);
      EXPECT_NEAR(energy_av(0, mc), -0.000605740233333333, 1e-8);
      EXPECT_NEAR(energy_av(1, mc), -0.030574223333333334, 0.03);
      EXPECT_NEAR(energy_av(2, mc), -0.089928316, 0.05);
      EXPECT_NEAR(energy_av(3, mc), -0.1784570533333333, 0.06);
      EXPECT_NEAR(energy_av(4, mc), -0.29619201333333334, 0.14);
      EXPECT_LE(mc.system().configuration().num_particles(), 5);

  //    // see if changing the c00 and c2N elements of colmat change the lnpi
  //    std::stringstream ss;
  //    criteria->bias().serialize(ss);
  //    TransitionMatrix tm(ss);
  //    tm.infrequent_update();
  //    INFO(feasst_str(tm.ln_prob().values()));
  ////    INFO("c00 " << tm.get_collection()->matrix()[0][0])
  ////    tm.get_collection()->increment(0, 0, 100);
  ////    INFO("c00 " << tm.get_collection()->matrix()[0][0])
  //    INFO("cn2 " << tm.get_collection()->matrix()[4][2])
  //    tm.get_collection()->increment(4, 2, 100);
  //    INFO("cn2 " << tm.get_collection()->matrix()[4][2])
  //    tm.infrequent_update();
  //    INFO(feasst_str(tm.ln_prob().values()));
    }
  }
}

TEST(MonteCarlo, spce_fh_LONG) {
  const std::vector<std::shared_ptr<Bias> > biases = {
    // MakeWangLandau({{"min_flatness", "30"}}),
    MakeTransitionMatrix({{"min_sweeps", "20"}}),
    MakeWLTM({{"collect_flatness", "20"},
              {"min_flatness", "25"},
              {"min_sweeps", "20"}})};
  for (auto bias : biases) {
    MonteCarlo mc;
    mc.set(spce({{"physical_constants", "CODATA2010"},
                 {"cubic_box_length", "20"},
                 {"alphaL", "5.6"},
                 {"kmax_squared", "38"}}));
    mc.set(MakeMetropolis({
      {"beta", str(1/kelvin2kJpermol(525, mc.configuration()))},
      {"chemical_potential", "-35.294567543492"}}));
    mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "0.275"}}));
    mc.add(MakeTrialRotate({{"weight", "1."}, {"tunable_param", "50."}}));
    mc.add(MakeTrialTransfer({{"particle_type", "0"}, {"weight", "4"}}));
    const double R = mc.configuration().physical_constants().ideal_gas_constant();
    const double temperature = 525*R/1000; // KJ/mol
    const double chemical_potential = -35.294567543492;
    auto criteria = MakeFlatHistogram(
      MakeMacrostateNumParticles(
        Histogram({{"width", "1"}, {"max", "5"}, {"min", "0"}})),
      bias,
      // MakeWangLandau({{"min_flatness", "30"}}),
      //MakeTransitionMatrix({{"min_sweeps", "20"}}),
      //MakeWLTM({{"collect_flatness", "15"},
      //          {"min_flatness", "20"},
      //          {"min_sweeps", "10"}})};
      {{"beta", str(1/temperature)},
       {"chemical_potential", str(chemical_potential)}});
    mc.set(criteria);
    const int steps_per = 1e5;
    mc.add(MakeLogAndMovie({{"steps_per", str(steps_per)}, {"file_name", "tmp/spce_fh"}}));
    mc.add(MakeCheckEnergyAndTune({{"steps_per", str(steps_per)}}));
    mc.add(MakeCriteriaUpdater({{"steps_per", str(steps_per)}}));
    mc.add(MakeCriteriaWriter({
      {"steps_per", str(steps_per)},
      {"file_name", "tmp/spce_crit.txt"},
    }));
    auto energy = MakeEnergy({
      {"file_name", "tmp/spce_fh_energy"},
      {"steps_per_update", "1"},
      {"steps_per_write", str(steps_per)},
      {"multistate", "true"}});
    mc.add(energy);
    mc.run_until_complete();
    test_serialize(mc);

    // known values of lnpi and energy
    const std::vector<std::vector<double> > lnpi_srsw = {
      {-2.8646724680467574586, 0.0131699993},
      {-1.9398627751910603179, 0.0117893901},
      {-1.5807236756003240075, 0.0163518954},
      {-1.5052807379855992487, 0.0066184482},
      {-1.5966420919511339349, 0.0104786041},
      {-1.7783427365938460074, 0.0291511326}};
    const std::vector<std::vector<double> >  en_srsw = {
      {0, 1e-13},
      {-0.0879115, 1.1293158298007674394e-05},
      {-2.25995, 0.038263},
      {-6.52141, 0.0519987},
      {-12.9855, 0.2504985},
      {-21.5192, 0.4465020}};

    INFO(bias->class_name());
    const LnProbability& lnpi = bias->ln_prob();
    for (int macro = 0; macro < lnpi.size(); ++macro) {
      EXPECT_NEAR(lnpi.value(macro), lnpi_srsw[macro][0],
        5*lnpi_srsw[macro][1]);
//      if (bias->class_name() == "TransitionMatrix") {
        const double en_std = std::sqrt(std::pow(en_srsw[macro][1], 2) +
          std::pow(energy->energy().block_stdev(), 2));
        EXPECT_NEAR(energy_av(macro, mc), en_srsw[macro][0], 10.*en_std);
//      }
    }

    EXPECT_LE(mc.system().configuration().num_particles(), 5);
  }
}

TEST(MonteCarlo, rpm_fh_LONG) {
  MonteCarlo mc;
  // mc.set(MakeRandomMT19937({{"seed", "default"}}));
  mc.set(MakeRandomMT19937({{"seed", "time"}}));
  mc.set(rpm({
    {"cubic_box_length", "12"},
    {"cutoff", "4.891304347826090"},
    {"alphaL", "6.87098396396261"}}));
  mc.add_to_reference(Potential(MakeDontVisitModel()));
  INFO("charge conversion " << CODATA2018().charge_conversion());
  const double temperature = 0.047899460618081;
  const double beta_mu = -13.94;
  auto criteria = MakeFlatHistogram(
    MakeMacrostateNumParticles(
      Histogram({{"width", "1"}, {"max", "2"}, {"min", "0"}}),
      {{"particle_type", "0"}}),
    // MakeWangLandau({{"min_flatness", "15"}}),
    MakeTransitionMatrix({{"min_sweeps", "100"}}),
    {{"beta", str(1/temperature)},
     {"chemical_potential0", str(beta_mu*temperature)},
     {"chemical_potential1", str(beta_mu*temperature)}});
  mc.set(criteria);
  INFO("beta_mu " << criteria->beta_mu(0));
  mc.add(MakeTrialTranslate({{"weight", "0.25"}, {"tunable_param", "0.1"}}));
  mc.add(MakeTrialTransferMultiple({
    {"weight", "4."},
    {"particle_type0", "0"},
    {"particle_type1", "1"},
    {"reference_index", "0"}}));
  const int steps_per = 1e5;
  if (criteria->bias().class_name() == "TransitionMatrix") {
    mc.add(MakeCriteriaUpdater({{"steps_per", str(steps_per)}}));
  }
  mc.add(MakeCriteriaWriter({
    {"steps_per", str(steps_per)},
    {"file_name", "tmp/rpmcrit.txt"}}));
  mc.add(MakeCheckProperties({{"steps_per", str(steps_per)}}));
  mc.add(MakeCheckPhysicality({{"steps_per", str(steps_per)}}));
  // mc.add(MakeCPUTime({{"steps_per", str(5*steps_per)}}));
  mc.add(MakeCheckNetCharge());
  mc.add(MakeLogAndMovie({{"steps_per", str(steps_per)}, {"file_name", "tmp/rpm_fh"}}));
  mc.add(MakeCheckEnergyAndTune({{"steps_per", str(steps_per)}}));
  mc.add(MakeEnergy({
    {"file_name", "tmp/rpm_fh_energy"},
    {"steps_per_update", "1"},
    {"steps_per_write", str(steps_per)},
    {"multistate", "true"}}));
  mc.run_until_complete();

  test_serialize(mc);

  const LnProbability& lnpi = criteria->bias().ln_prob();
  EXPECT_NEAR(lnpi.value(0), -1.2994315780357, 0.1);
  EXPECT_NEAR(lnpi.value(1), -1.08646312498868, 0.1);
  EXPECT_NEAR(lnpi.value(2), -0.941850889679828, 0.1);
  EXPECT_NEAR(energy_av(0, mc), 0, 1e-14);
  EXPECT_NEAR(energy_av(1, mc), -0.939408, 0.02);
  EXPECT_NEAR(energy_av(2, mc), -2.02625, 0.04);
}

TEST(MonteCarlo, rpm_fh_divalent_LONG) {
  const double temperature = 0.25;
  //const double temperature = 0.05;
  const double beta_mu = -7.94;
  //const double beta_mu = -23.94;
  MonteCarlo mc;
  mc.set(MakeRandomMT19937({{"seed", "time"}}));
  mc.set(rpm({
    {"delta", "0.3"},
    {"charge_ratio", "2"},
    {"cubic_box_length", "15"},
    {"cutoff", "7.5"},
    {"kmax_squared", "25"},
    {"alphaL", "5"}}));
  mc.add_to_reference(Potential(MakeDontVisitModel()));
//  mc.set(1, Potential(MakeModelTwoBodyFactory({MakeHardSphere(),
//                                               MakeChargeScreened()}),
//                      MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
  auto criteria = MakeFlatHistogram(
    MakeMacrostateNumParticles(
      Histogram({{"width", "1"}, {"max", "5"}, {"min", "0"}}),
      {{"particle_type", "0"}}),
    // MakeWangLandau({{"min_flatness", "100"}}),
    MakeTransitionMatrix({{"min_sweeps", "100"}}),
    {{"beta", str(1/temperature)},
     {"chemical_potential0", str(beta_mu*temperature)},
     {"chemical_potential1", str(beta_mu*temperature)}});
  mc.set(criteria);
  mc.add(MakeTrialTranslate({{"weight", "0.25"}, {"tunable_param", "0.1"}}));
  mc.add(MakeTrialTransferMultiple({
    {"weight", "1."},
    {"particle_type0", "0"},
    {"particle_type1", "1"},
    {"particle_type2", "1"},
    {"reference_index", "0"}}));
//  auto neighbor_criteria = MakeNeighborCriteria({{"maximum_distance", "1.5"},
//                                                 {"minimum_distance", "1"},
//                                                 {"site_type0", "0"},
//                                                 {"site_type1", "1"},
//                                                 {"potential_index", "1"}});
//  add_rigid_cluster_trials(&mc,
//    neighbor_criteria,
//    {{"tunable_param", "50"}});
  const int steps_per = 1e5;
  mc.add(MakeCriteriaUpdater({{"steps_per", str(steps_per)}}));
  mc.add(MakeCriteriaWriter({
    {"steps_per", str(steps_per)},
    {"file_name", "tmp/dival_fh_crit.txt"}}));
  mc.add(MakeLogAndMovie({{"steps_per", str(steps_per)}, {"file_name", "tmp/dival_fh"}}));
  mc.add(MakeCheckEnergyAndTune({{"steps_per", str(steps_per)}, {"tolerance", str(1e-4)}}));
  mc.add(MakeCheckNetCharge({{"steps_per", str(steps_per)}}));
  const int en_index = mc.num_analyzers();
  mc.add(MakeEnergy({
    {"file_name", "tmp/dival_fh_energy"},
    {"steps_per_update", "1"},
    {"steps_per_write", str(steps_per)},
    {"multistate", "true"}}));
  mc.run_until_complete();
  // mc.attempt(1e7);

  const LnProbability lnpi = FlatHistogram(mc.criteria()).bias().ln_prob();
  EXPECT_NEAR(lnpi.value(0), -6.6615, 0.06);
  EXPECT_NEAR(lnpi.value(1), -3.6256, 0.06);
  EXPECT_NEAR(lnpi.value(2), -2.1046, 0.06);
  EXPECT_NEAR(lnpi.value(3), -1.3685, 0.06);
  EXPECT_NEAR(lnpi.value(4), -1.1371, 0.06);
  EXPECT_NEAR(lnpi.value(5), -1.2911, 0.06);
  const std::vector<std::shared_ptr<Analyze> >& en =
    mc.analyzers()[en_index]->analyzers();
  EXPECT_NEAR(en[0]->accumulator().average(), 0, 1e-14);
  EXPECT_NEAR(en[1]->accumulator().average(), -1.30701, 0.03);
  EXPECT_NEAR(en[2]->accumulator().average(), -2.98115, 0.03);
  EXPECT_NEAR(en[3]->accumulator().average(), -4.85254, 0.05);
  EXPECT_NEAR(en[4]->accumulator().average(), -6.80956, 0.12);
  EXPECT_NEAR(en[5]->accumulator().average(), -8.85025, 0.16);
}

}  // namespace feasst
