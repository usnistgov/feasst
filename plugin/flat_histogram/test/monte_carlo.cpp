#include "utils/test/utils.h"
#include "math/include/histogram.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/domain.h"
#include "system/include/hard_sphere.h"
#include "system/include/dont_visit_model.h"
#include "system/include/model_two_body_factory.h"
#include "system/include/potential.h"
#include "system/include/lennard_jones.h"
#include "system/include/long_range_corrections.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/run.h"
#include "monte_carlo/include/trial_transfer.h"
#include "monte_carlo/include/trial_rotate.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/trial_add.h"
#include "steppers/include/energy.h"
#include "steppers/include/criteria_writer.h"
#include "steppers/include/criteria_updater.h"
#include "steppers/include/tune.h"
#include "steppers/include/cpu_time.h"
#include "steppers/include/check_properties.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/check_physicality.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/tune.h"
#include "steppers/include/log_and_movie.h"
#include "charge/include/ewald.h"
#include "charge/include/charge_screened.h"
#include "charge/include/charge_self.h"
#include "charge/include/check_net_charge.h"
#include "charge/include/utils.h"
#include "charge/test/charge_utils.h"
#include "charge/include/trial_transfer_multiple.h"
#include "flat_histogram/include/flat_histogram.h"
#include "flat_histogram/include/macrostate_num_particles.h"
#include "flat_histogram/include/transition_matrix.h"
#include "flat_histogram/include/wang_landau.h"
#include "flat_histogram/include/wltm.h"
#include "flat_histogram/include/ensemble.h"

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
TEST(MonteCarlo, ideal_gas_fh_eos_LONG) {
  MonteCarlo monte_carlo;
  monte_carlo.add(MakeConfiguration({{"cubic_box_length", "8"},
    {"particle_type", install_dir() + "/forcefield/atom.fstprt"}}));
  monte_carlo.add(MakePotential(MakeDontVisitModel()));
  monte_carlo.set(MakeThermoParams({{"beta", str(1./1.2)}, {"chemical_potential", "-3"}}));
//  auto criteria = MakeFlatHistogram(
//      MakeMacrostateNumParticles(Histogram({{"width", "1"}, {"min", "0"}, {"max", "50"}})),
//      MakeTransitionMatrix({{"min_sweeps", "100"}}));
//  auto criteria = MakeFlatHistogram(
//      MakeMacrostateNumParticles({{"width", "1"}, {"min", "0"}, {"max", "50"}}),
//      MakeTransitionMatrix({{"min_sweeps", "100"}}));
  auto criteria = MakeFlatHistogram({
    {"Macrostate", "MacrostateNumParticles"}, {"width", "1"}, {"min", "0"}, {"max", "50"},
    {"Bias", "TransitionMatrix"}, {"min_sweeps", "100"}});
  monte_carlo.set(criteria);
  monte_carlo.add(MakeTrialTransfer({{"particle_type", "0"}}));
  monte_carlo.add(MakeCriteriaUpdater({{"trials_per", str(1e5)}}));
  monte_carlo.add(MakeCriteriaWriter({{"trials_per", str(1e5)}, {"file_name", "tmp/id_fh.txt"}}));
  monte_carlo.run_until_complete();
  GrandCanonicalEnsemble gce(*criteria, monte_carlo.system().thermo_params().beta_mu());
  for (double delta_conjugate = -6; delta_conjugate < 1; delta_conjugate += 0.1) {
    gce.reweight(delta_conjugate);
    if (gce.ln_prob().value(gce.ln_prob().size() - 1) < -6) {
      EXPECT_NEAR(gce.average_macrostate()/gce.betaPV(), 1, 0.01);
    }
  }
}

TEST(MonteCarlo, hard_sphere_LONG) {
  MonteCarlo mc;
  //mc.set(MakeRandomMT19937({{"seed", "123"}}));
  mc.add(MakeConfiguration({{"cubic_box_length", "8"},
    {"particle_type", "../forcefield/hard_sphere.fstprt"}}));
  mc.add(MakePotential(MakeHardSphere()));
//  mc.add_to_optimized(MakePotential(MakeHardSphere(), MakeVisitModelCell({{"min_length", "1"}})));
  mc.set(MakeThermoParams({{"beta", "1"}, {"chemical_potential", "-2.352321"}}));
  //auto bias = MakeWLTM({{"collect_flatness", "18"},
  //              {"min_flatness", "22"},
  //              {"min_sweeps", "100"}});
  auto bias = MakeTransitionMatrix({{"min_sweeps", "100"}});
  mc.set(MakeFlatHistogram(
      MakeMacrostateNumParticles(
          Histogram({{"width", "1"}, {"max", "40"}, {"min", "0"}})), bias));
  const std::string is_new_only = "true";
  mc.add(MakeTrialTranslate({{"new_only", is_new_only}, {"weight", "1."}, {"tunable_param", "1."}}));
  mc.add(MakeTrialTransfer({{"weight", "4"}, {"particle_type", "0"}}));
  const std::string trials_per = "100";
  mc.add(MakeCheckEnergy({{"trials_per", trials_per}, {"tolerance", "0.0001"}}));
  mc.add(MakeCheckPhysicality({{"trials_per", "1"}}));
  mc.add(MakeTune({{"stop_after_phase", "0"}}));
  mc.add(MakeLogAndMovie({{"trials_per", trials_per},
                          {"file_name", "hs_fh"},
                          {"file_name_append_phase", "True"}}));
  mc.add(MakeCriteriaUpdater({{"trials_per", trials_per}}));
  mc.add(MakeCriteriaWriter({{"trials_per", trials_per},
                             {"file_name", "tmp/crit.txt"},
                             {"file_name_append_phase", "True"}}));
  mc.run_until_complete();
  INFO(feasst_str(bias->ln_prob().values()));
  EXPECT_NEAR(bias->ln_prob().value(0), -41.16903361558974, 0.25);
  //EXPECT_NEAR(bias->ln_prob().value(0), -41.3327752, 0.05);
}

MonteCarlo test_lj_fh(const int num_steps,
    const std::string bias_name,
    int sweeps = 10,
    bool test_multi = false,
    const int min = 1,
    const int max = 5,
    bool dont_use_multi = false) {
  MonteCarlo mc;
  int ref = -1;
  mc.add(MakeConfiguration({{"cubic_box_length", "8"},
                            {"particle_type0", "../forcefield/lj.fstprt"}}));
  mc.add(MakePotential(MakeLennardJones()));
  mc.add(MakePotential(MakeLongRangeCorrections()));
  if (num_steps == 1) {
    if (test_multi) {
      ref = 0;
      mc.add_to_reference(MakePotential(MakeDontVisitModel()));
    }
  } else {
    ref = 0;
    mc.run(MakeConvertToRefPotential({{"cutoff", "1"}, {"use_cell", "true"}}));
  }
  mc.set(MakeThermoParams({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc.add(MakeTrialAdd({{"particle_type", "0"}}));
  mc.run(MakeRun({{"until_num_particles", str(min)}}));
  mc.run(MakeRemoveTrial({{"name", "TrialAdd"}}));
  argtype transfer_args =
    { {"particle_type0", "0"},
      {"reference_index", str(ref)},
      {"num_steps", str(num_steps)},
      {"weight", "4"}};
  if (test_multi) {
    transfer_args.insert({"particle_type1", "0"});
    transfer_args.insert({"shift", "-2"});
  }
  if (dont_use_multi) {
    ASSERT(!test_multi, "er");
    mc.add(MakeTrialTransfer({ {"particle_type", "0"},
      {"reference_index", str(ref)},
      {"num_steps", str(num_steps)},
      {"weight", "4"}}));
  } else {
    mc.add(MakeTrialTransferMultiple(transfer_args));
  }
  EXPECT_EQ(mc.trial(0).weight(), 1);
  EXPECT_EQ(mc.trial(1).weight(), 2);
  EXPECT_EQ(mc.trial(2).weight(), 2);
  std::shared_ptr<Bias> bias;
  if (bias_name == "TM") {
    bias = MakeTransitionMatrix({{"min_sweeps", str(sweeps)}});
    //bias = MakeTransitionMatrix({{"min_sweeps", str(1e5)}, {"new_sweep", "1"}});
  } else if (bias_name == "WL") {
    if (sweeps == 10) sweeps = 20;
    bias = MakeWangLandau({{"min_flatness", str(sweeps)}});
  } else if (bias_name == "WLTM") {
    bias = MakeWLTM({{"collect_flatness", "15"},
      {"min_flatness", "20"}, {"min_sweeps", str(sweeps)}});
  } else {
    FATAL("unrecognized");
  }
  std::string width = "1";
  if (test_multi) width = "2";
  mc.set(MakeThermoParams({{"beta", str(1./1.5)},
      {"chemical_potential", "-2.352321"}}));
//      {{"soft_max", "5"}, {"soft_min", "1"}}));
//      {{"particle_type", "0"}}));
  auto criteria = MakeFlatHistogram(
    MakeMacrostateNumParticles(
      Histogram({{"width", width}, {"max", str(max)}, {"min", str(min)}})),
    bias);
  INFO(criteria->bias().class_name());
  mc.set(criteria);
  const std::string trials_per(str(1e3));
  // const std::string trials_per(str(1e4));
  mc.add(MakeLogAndMovie({{"trials_per", str(trials_per)}, {"file_name", "tmp/lj_fh"}}));
  mc.add(MakeCheckEnergy({{"trials_per", str(trials_per)}}));
  //mc.add(MakeCheckEnergyAndTune({{"trials_per", str(trials_per)}}));
  mc.add(MakeTune({{"multistate", "true"}, {"trials_per_write", str(trials_per)}, {"file_name", "tmp/tune.txt"}}));
  mc.add(MakeCriteriaUpdater({{"trials_per", str(1)}}));
  mc.add(MakeCriteriaWriter({
    {"trials_per", trials_per},
    {"file_name", "tmp/ljcrit.txt"},
    {"file_name_append_phase", "true"}}));
  auto energy = MakeEnergy({
    {"file_name", "tmp/lj_fh_energy"},
    {"file_name_append_phase", "true"},
    {"trials_per_update", "1"},
    {"trials_per_write", str(trials_per)},
    {"multistate", "true"}});
  EXPECT_EQ(energy->trials_per_update(), 1);
  EXPECT_EQ(energy->trials_per_write(), 1e3);
  mc.add(energy);
  return mc;
}

TEST(MonteCarlo, lj_fh_01) {
  MonteCarlo mc = test_lj_fh(1, "TM", 10, false, 0, 1);
  mc.run_until_complete();
  FlatHistogram fh(mc.criteria());
  const LnProbability lnpi = fh.bias().ln_prob();
  EXPECT_NEAR(lnpi.value(1) - lnpi.value(0), 4.67, 0.2);

  // obtain tm/cm
  std::stringstream ss;
  fh.bias().serialize(ss);
  TransitionMatrix tm(ss);
  INFO(tm.collection().min_blocks());
  std::vector<LnProbability> ln_probs = tm.collection().ln_prob_blocks();
  Accumulator acc;
  for (const auto& ln_prob : ln_probs) {
    acc.accumulate(ln_prob.value(1) - ln_prob.value(0));
//    INFO(ln_prob.value(1) - ln_prob.value(0));
    //INFO(feasst_str(ln_prob.values()));
  }
//  INFO(acc.stdev_of_av());
//  INFO(fh.write());
}

TEST(MonteCarlo, lj_fh_10sweep_LONG) {
  //for (int num_steps : {1}) {
  //for (int num_steps : {2}) {
  for (int num_steps : {1, 2}) {
    //for (const std::string bias_name : {"TM"}) {
    for (const std::string bias_name : {"TM", "WL", "WLTM"}) {
      MonteCarlo mc = test_serialize(test_lj_fh(num_steps, bias_name));
      //mc.attempt(1e4);
      //mc.attempt(1e5); // note more than 1e4 steps required for TM
      mc.run_until_complete();
      // INFO(mc.criteria().write());

      // compare with known values of lnpi
      //const LnProbability * lnpi = &criteria->bias().ln_prob();
      const LnProbability lnpi = FlatHistogram(mc.criteria()).bias().ln_prob();
      //EXPECT_NEAR(lnpi.value(0), -18.707570324988800000, 0.55);
      EXPECT_NEAR(lnpi.value(0), -14.037373358321800000, 0.75);
      EXPECT_NEAR(lnpi.value(1), -10.050312091655200000, 0.6);
      EXPECT_NEAR(lnpi.value(2), -6.458920624988570000, 0.55);
      EXPECT_NEAR(lnpi.value(3), -3.145637424988510000, 0.55);
      EXPECT_NEAR(lnpi.value(4), -0.045677458321876000, 0.55);

      // compare with known values of energy
      //EXPECT_NEAR(energy_av(0, mc), 0, 1e-14);
      EXPECT_NEAR(energy_av(0, mc), -0.000605740233333333, 1e-8);
      EXPECT_NEAR(energy_av(1, mc), -0.030574223333333334, 0.05);
      EXPECT_NEAR(energy_av(2, mc), -0.089928316, 0.08);
      EXPECT_NEAR(energy_av(3, mc), -0.1784570533333333, 0.11);
      EXPECT_NEAR(energy_av(4, mc), -0.29619201333333334, 0.15);
      EXPECT_LE(mc.system().configuration().num_particles(), 5);
    }
  }
}

TEST(MonteCarlo, lj_fh_block) {
  MonteCarlo mc = test_serialize(test_lj_fh(1, "TM", 10));
  mc.run_until_complete();
//  const LnProbability lnpi = FlatHistogram(mc.criteria()).bias().ln_prob();
}

TEST(MonteCarlo, soft_min_macro) {
  auto mc = MakeMonteCarlo({{
    {"Configuration", {{"particle_type0", "../forcefield/lj.fstprt"},
                       {"cubic_box_length", "8"}}},
    {"Potential", {{"Model", "LennardJones"}}},
    {"ThermoParams", {{"beta", "1.2"}, {"chemical_potential", "1."}}},
    {"FlatHistogram", {{"Macrostate", "MacrostateNumParticles"}, {"width", "1"},
      {"max", "5"}, {"min", "0"}, {"soft_macro_min", "1"}, {"soft_macro_max", "4"},
      {"Bias", "TransitionMatrix"}, {"min_sweeps", "10"}}},
  }});
  const FlatHistogram& fh = FlatHistogram(mc->criteria());
  EXPECT_NEAR(1, fh.macrostate().value(0), NEAR_ZERO);
}

TEST(MonteCarlo, lj_fh_with0) {
  //for (int num_steps : {1}) {
  for (int num_steps : {1, 2}) {
    //for (const std::string bias_name : {"WLTM"}) {
    for (const std::string bias_name : {"TM", "WL", "WLTM"}) {
      for (const bool dont_use_multi : {true, false}) {
        MonteCarlo mc = test_serialize(test_lj_fh(num_steps, bias_name, 10, false, 0, 5, dont_use_multi));
        //mc.attempt(1e4);
        //mc.attempt(1e5); // note more than 1e4 steps required for TM
        mc.run_until_complete();
        // INFO(mc.criteria().write());

        // compare with known values of lnpi
        //const LnProbability * lnpi = &criteria->bias().ln_prob();
        const LnProbability lnpi = FlatHistogram(mc.criteria()).bias().ln_prob();
        //INFO(feasst_str(lnpi.values()));
        EXPECT_NEAR(lnpi.value(0), -18.707570324988800000, 0.55);
        EXPECT_NEAR(lnpi.value(1), -14.037373358321800000, 0.75);
        EXPECT_NEAR(lnpi.value(2), -10.050312091655200000, 0.6);
        EXPECT_NEAR(lnpi.value(3), -6.458920624988570000, 0.55);
        EXPECT_NEAR(lnpi.value(4), -3.145637424988510000, 0.55);
        EXPECT_NEAR(lnpi.value(5), -0.045677458321876000, 0.55);

        // compare with known values of energy
        EXPECT_NEAR(energy_av(0, mc), 0, 1e-14);
        EXPECT_NEAR(energy_av(1, mc), -0.000605740233333333, 1e-8);
        EXPECT_NEAR(energy_av(2, mc), -0.030574223333333334, 0.05);
        EXPECT_NEAR(energy_av(3, mc), -0.089928316, 0.08);
        EXPECT_NEAR(energy_av(4, mc), -0.1784570533333333, 0.11);
        EXPECT_NEAR(energy_av(5, mc), -0.29619201333333334, 0.15);
        EXPECT_LE(mc.system().configuration().num_particles(), 5);
      }
    }
  }
}

TEST(MonteCarlo, lj_fh_LONG) {
  MonteCarlo mc = test_serialize(test_lj_fh(4, "TM", 1000));
  mc.run_until_complete();
  const LnProbability lnpi = FlatHistogram(mc.criteria()).bias().ln_prob();
  EXPECT_NEAR(lnpi.value(0), -14.037373358321800000, 0.02);
  EXPECT_NEAR(lnpi.value(1), -10.050312091655200000, 0.02);
  EXPECT_NEAR(lnpi.value(2), -6.458920624988570000, 0.02);
  EXPECT_NEAR(lnpi.value(3), -3.145637424988510000, 0.01);
  EXPECT_NEAR(lnpi.value(4), -0.045677458321876000, 0.005);
  EXPECT_NEAR(energy_av(0, mc), -0.000605740233333333, 1e-8);
  EXPECT_NEAR(energy_av(1, mc), -0.030574223333333334, 0.001);
  EXPECT_NEAR(energy_av(2, mc), -0.089928316, 0.003);
  EXPECT_NEAR(energy_av(3, mc), -0.1784570533333333, 0.005);
  EXPECT_NEAR(energy_av(4, mc), -0.29619201333333334, 0.0075);
  const LnProbability lnpi3 = lnpi.reduce(2);
  INFO(feasst_str(lnpi3.values()));
  EXPECT_NEAR(lnpi3.value(0), -13.9933350923078, 0.025);
  EXPECT_NEAR(lnpi3.value(1), -6.41488235897456, 0.025);
  EXPECT_NEAR(lnpi3.value(2), -0.00163919230786818, 0.005);
}

// HWH Test for fixing dccb translate issue
//TEST(MonteCarlo, lj_fh_transition_VERY_LONG) {
//  for (int num_steps : {1, 4}) {
//    MonteCarlo mc = test_serialize(test_lj_fh(num_steps, "TM", 10000, false, 150, 151));
//    mc.set(MakeThermoParams({{"beta", str(1./0.75)}, {"chemical_potential", "-4.05045075"}}));
//    mc.run_until_complete();
//    const LnProbability lnpi = FlatHistogram(mc.criteria()).bias().ln_prob();
//    INFO("num_steps: " << num_steps << " " << feasst_str(lnpi.values()));
//    EXPECT_NEAR(lnpi.value(1) - lnpi.value(0), -241.79257 - -241.79008, 0.005);
//    //EXPECT_NEAR(energy_av(0, mc), 0, 0.5);
//    //EXPECT_NEAR(energy_av(1, mc), 0, 0.5);
//  }
//}

TEST(MonteCarlo, lj_fh_liquid_LONG) {
  MonteCarlo mc = test_serialize(test_lj_fh(4, "TM", 1000, false, 100, 105));
  mc.run_until_complete();
  const LnProbability lnpi = FlatHistogram(mc.criteria()).bias().ln_prob();
  EXPECT_NEAR(lnpi.value(0), -4.92194963175925, 0.025);
  EXPECT_NEAR(lnpi.value(1), -4.03855513175926, 0.02);
  EXPECT_NEAR(lnpi.value(2), -3.15822813175925, 0.02);
  EXPECT_NEAR(lnpi.value(3), -2.28019483175925, 0.015);
  EXPECT_NEAR(lnpi.value(4), -1.40647303175926, 0.0075);
  EXPECT_NEAR(lnpi.value(5), -0.535594831759248, 0.005);
  EXPECT_NEAR(energy_av(0, mc), -1.381223800E+02, 0.5);
  EXPECT_NEAR(energy_av(1, mc), -1.408257000E+02, 0.5);
  EXPECT_NEAR(energy_av(2, mc), -1.435426000E+02, 0.5);
  EXPECT_NEAR(energy_av(3, mc), -1.462802100E+02, 0.5);
  EXPECT_NEAR(energy_av(4, mc), -1.490501400E+02, 0.5);
  EXPECT_NEAR(energy_av(5, mc), -1.518517300E+02, 0.6);
}

TEST(MonteCarlo, lj_fh_multi_LONG) {
  //MonteCarlo mc = test_serialize(test_lj_fh(1, "WLTM", 1000, true));
  //MonteCarlo mc = test_serialize(test_lj_fh(1, "WL", 200, true));
  MonteCarlo mc = test_serialize(test_lj_fh(1, "TM", 1000, true));
  //mc.set(MakeRandomMT19937({{"seed", "123"}}));
  mc.run_until_complete();
  INFO(mc.criteria().write());
  const LnProbability lnpi = FlatHistogram(mc.criteria()).bias().ln_prob();
  EXPECT_NEAR(lnpi.value(0), -13.9933350923078, 0.0225);
  EXPECT_NEAR(lnpi.value(1), -6.41488235897456, 0.02);
  EXPECT_NEAR(lnpi.value(2), -0.00163919230786818, 0.005);
  EXPECT_NEAR(energy_av(0, mc), -0.000605740233333333, 1e-8);
  EXPECT_NEAR(energy_av(1, mc), -0.089928316, 0.0025);
  EXPECT_NEAR(energy_av(2, mc), -0.29619201333333334, 0.02);
}

// HWH add num steps to spce fh test for DCCB diagnosis
MonteCarlo test_spce_fh(std::shared_ptr<Bias> bias,
    const int num_steps = 1,
    bool test = true,
    const int min = 0,
    const int max = 5,
    const int trials_per = 1e3) {
  INFO(bias->class_name());
  MonteCarlo mc;
  argtype spce_args = {{"physical_constants", "CODATA2010"},
                       {"cubic_box_length", "20"},
                       {"alpha", str(5.6/20)},
                       {"kmax_squared", "38"}};
  int ref = -1;
  if (num_steps > 1) {
    //spce_args.insert({"dual_cut", str(10)});
    spce_args.insert({"dual_cut", str(3.16555789)});
    ref = 0;
  }
  mc.set(spce(spce_args));
  const double beta = 1/kelvin2kJpermol(525, mc.configuration()); // mol/kJ
  mc.set(MakeThermoParams({{"beta", "1"}, {"chemical_potential", "1"}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "0.275"}}));
  mc.add(MakeTrialRotate({{"weight", "1."}, {"tunable_param", "0.2"}}));
  mc.add(MakeTrialAdd({{"particle_type", "0"}}));
  mc.run(MakeRun({{"until_num_particles", str(min)}}));
  mc.run(MakeRemoveTrial({{"name", "TrialAdd"}}));
  mc.set(MakeThermoParams({{"beta", str(beta)}, {"chemical_potential", str(-8.14/beta)}}));
  mc.add(MakeTrialTransfer({
    {"particle_type", "0"},
    {"weight", "4"},
    {"reference_index", str(ref)},
    {"num_steps", str(num_steps)}}));
  auto criteria = MakeFlatHistogram(
    MakeMacrostateNumParticles(
      Histogram({{"width", "1"}, {"max", str(max)}, {"min", str(min)}})),
    bias);
  mc.set(criteria);
  mc.add(MakeLogAndMovie({{"trials_per", str(trials_per)}, {"file_name", "tmp/spce_fh"}}));
  mc.add(MakeCheckEnergy({{"trials_per", str(trials_per)}, {"tolerance", str(1e-6)}}));
  mc.add(MakeTune({{"multistate", "true"}, {"trials_per_write", str(trials_per)}, {"file_name", "tmp/spce_tune.txt"}}));
  mc.add(MakeCriteriaUpdater({{"trials_per", str(trials_per)}}));
  mc.add(MakeCriteriaWriter({
    {"trials_per", str(trials_per)},
    {"file_name", "tmp/spce_crit.txt"}}));
  auto energy = MakeEnergy({
    {"file_name", "tmp/spce_fh_energy"},
    {"trials_per_update", "1"},
    {"trials_per_write", str(trials_per)},
    {"multistate", "true"}});
  mc.add(energy);
  MonteCarlo mc2 = test_serialize(mc);
  mc2.run_until_complete();

  if (!test) return mc2;

  EXPECT_LE(mc2.system().configuration().num_particles(), 5);

  // known values of lnpi and energy
  const std::vector<std::vector<double> > lnpi_srsw = {
    {-2.7207, 0.015},
    {-1.8523, 0.015},
    {-1.54708, 0.016},
    {-1.51786, 0.015},
    {-1.6479, 0.015},
    {-1.8786, 0.03}};
  const std::vector<std::vector<double> >  en_srsw = {
    {0, 1e-13},
    {-0.0879115, 1.1293158298007674394e-05},
    {-2.326, 0.12},
    {-6.806, 0.24},
    {-13.499, 0.5},
    {-22.27, 1.0}};

  FlatHistogram fh(mc2.criteria());
  const LnProbability& lnpi = fh.bias().ln_prob();
  for (int macro = 0; macro < lnpi.size(); ++macro) {
    EXPECT_NEAR(lnpi.value(macro), lnpi_srsw[macro][0],
      15*lnpi_srsw[macro][1]);
//      if (bias->class_name() == "TransitionMatrix") {
      const double en_std = std::sqrt(std::pow(en_srsw[macro][1], 2) +
        std::pow(energy->energy().block_stdev(), 2));
      EXPECT_NEAR(energy_av(macro, mc2), en_srsw[macro][0], 15.*en_std);
//      }
  }

  return mc2;
}

TEST(MonteCarlo, spce_fh_LONG) {
  const std::vector<std::shared_ptr<Bias> > biases = {
    MakeWangLandau({{"min_flatness", "40"}}),
    MakeTransitionMatrix({{"min_sweeps", "25"}}),
    MakeWLTM({{"collect_flatness", "20"},
              {"min_flatness", "25"},
              {"min_sweeps", "20"}})};
  for (auto bias : biases) test_spce_fh(bias);
}

TEST(MonteCarlo, spce_fh_VERY_LONG) {
  //for (int num_steps : {1}) {
  for (int num_steps : {4}) {
  //for (int num_steps : {1, 4}) {
    MonteCarlo mc = test_spce_fh(
      MakeTransitionMatrix({{"min_sweeps", "1000"}}),
      num_steps,
      false); // test
    FlatHistogram fh(mc.criteria());
    INFO(feasst_str(fh.bias().ln_prob().values()));
    const LnProbability& lnpi = fh.bias().ln_prob();
    EXPECT_NEAR(lnpi.value(0), -2.72070275309203, 0.02);
    EXPECT_NEAR(lnpi.value(1), -1.85234049431879, 0.02);
    EXPECT_NEAR(lnpi.value(2), -1.54708325224374, 0.02);
    EXPECT_NEAR(lnpi.value(3), -1.51786213939762, 0.02);
    EXPECT_NEAR(lnpi.value(4), -1.64791755404893, 0.02);
    EXPECT_NEAR(lnpi.value(5), -1.87860075480337, 0.05);
    const std::vector<std::shared_ptr<Analyze> >& en =
      mc.analyzers().back()->analyzers();
    EXPECT_NEAR(en[0]->accumulator().average(), 0, 1e-13);
    EXPECT_NEAR(en[1]->accumulator().average(), -0.08790895, 1e-5);
    EXPECT_NEAR(en[2]->accumulator().average(), -2.32656326, 0.2);
    EXPECT_NEAR(en[3]->accumulator().average(), -6.80645068, 0.5);
    EXPECT_NEAR(en[4]->accumulator().average(), -13.49913788, 1.);
    EXPECT_NEAR(en[5]->accumulator().average(), -22.27407753, 2.25);
  }
}

// This one shows difference with num_steps. Update lnpi and en checks
TEST(MonteCarlo, spce_fh_liquid_VERY_LONG) {
  //for (int num_steps : {1}) {
  //for (int num_steps : {2}) {
  for (int num_steps : {4}) {
  //for (int num_steps : {16}) {
  //for (int num_steps : {1, 4}) {
    MonteCarlo mc = test_spce_fh(
      MakeTransitionMatrix({{"min_sweeps", "100"}}),
      num_steps,
      false, // test
      100,   // min
      105,   // max
      1e4);  // trials_per
    FlatHistogram fh(mc.criteria());
    INFO(feasst_str(fh.bias().ln_prob().values()));
    const LnProbability& lnpi = fh.bias().ln_prob();
    EXPECT_NEAR(lnpi.value(0), -1.9471154, 0.15);
    EXPECT_NEAR(lnpi.value(1), -1.898168, 0.15);
    EXPECT_NEAR(lnpi.value(2), -1.8426095, 0.15);
    EXPECT_NEAR(lnpi.value(3), -1.76756244, 0.15);
    EXPECT_NEAR(lnpi.value(4), -1.707793039, 0.15);
    EXPECT_NEAR(lnpi.value(5), -1.62427498, 0.15);
    const std::vector<std::shared_ptr<Analyze> >& en =
      mc.analyzers().back()->analyzers();
    EXPECT_NEAR(en[0]->accumulator().average(), -2666.66, 50);
    EXPECT_NEAR(en[1]->accumulator().average(), -2703.74, 50);
    EXPECT_NEAR(en[2]->accumulator().average(), -2738.90, 50);
    EXPECT_NEAR(en[3]->accumulator().average(), -2772.03, 50);
    EXPECT_NEAR(en[4]->accumulator().average(), -2804.02, 50);
    EXPECT_NEAR(en[5]->accumulator().average(), -2836.26, 50);
  }
}

MonteCarlo rpm_fh_test(
    const int min = 0,
    const int max = 2,
    const int trials_per = 1e3,
    const int num_steps = 1) {
  MonteCarlo mc;
  // mc.set(MakeRandomMT19937({{"seed", "default"}}));
  mc.set(MakeRandomMT19937({{"seed", "time"}}));
  argtype rpm_args = {{"cubic_box_length", "12"},
    {"cutoff", "4.891304347826090"},
    {"kmax_squared", "38"},
    {"alpha", str(6.87098396396261/12)}};
  int ref = -1;
  if (num_steps > 1) {
    rpm_args.insert({"dual_cut", "1"});
    ref = 0;
  }
  mc.add(MakeConfiguration({{"cubic_box_length", "12"},
    {"particle_type0", install_dir() + "/plugin/charge/forcefield/rpm_plus.fstprt"},
    {"particle_type1", install_dir() + "/plugin/charge/forcefield/rpm_minus.fstprt"},
    {"cutoff", "4.891304347826090"},
    {"charge0", str( 1/std::sqrt(CODATA2018().charge_conversion()))},
    {"charge1", str(-1/std::sqrt(CODATA2018().charge_conversion()))}}));
  mc.add(MakePotential(MakeEwald({{"alpha", str(6.87098396396261/12)},
    {"kmax_squared", "38"}})));
  mc.add(MakePotential(MakeModelTwoBodyFactory(MakeHardSphere(), MakeChargeScreened()),
                                    {{"table_size", "1e6"}}));
  mc.add(MakePotential(MakeChargeSelf()));

  mc.add_to_reference(MakePotential(MakeDontVisitModel()));
  INFO("charge conversion " << MAX_PRECISION << CODATA2018().charge_conversion());
  const double temperature = 0.047899460618081;
  const double beta_mu = -13.94;
  mc.set(MakeThermoParams({{"beta", str(1/temperature)},
     {"chemical_potential0", str(beta_mu*temperature)},
     {"chemical_potential1", str(beta_mu*temperature)}}));
  auto criteria = MakeFlatHistogram(
    MakeMacrostateNumParticles(
      Histogram({{"width", "1"}, {"max", "2"}, {"min", "0"}}),
      {{"particle_type", "0"}}),
    // MakeWangLandau({{"min_flatness", "15"}}),
    MakeTransitionMatrix({{"min_sweeps", "100"}}));
  mc.set(criteria);
  INFO("beta_mu " << mc.system().thermo_params().beta_mu(0));
  mc.add(MakeTrialTranslate({
    {"weight", "0.25"},
    {"tunable_param", "0.1"},
    {"reference_index", str(ref)},
    {"num_steps", str(num_steps)}}));
  mc.add(MakeTrialTransferMultiple({
    {"weight", "4."},
    {"particle_type0", "0"},
    {"particle_type1", "1"},
    {"reference_index", "0"},
    {"num_steps", str(num_steps)}}));
  if (criteria->bias().class_name() == "TransitionMatrix") {
    mc.add(MakeCriteriaUpdater({{"trials_per", str(trials_per)}}));
  }
  mc.add(MakeCriteriaWriter({
    {"trials_per", str(trials_per)},
    {"file_name", "tmp/rpmcrit.txt"}}));
  mc.add(MakeCheckProperties({{"trials_per", str(trials_per)}}));
  mc.add(MakeCheckPhysicality({{"trials_per", str(trials_per)}}));
  // mc.add(MakeCPUTime({{"trials_per", str(5*trials_per)}}));
  mc.add(MakeCheckNetCharge());
  mc.add(MakeLogAndMovie({{"trials_per", str(trials_per)}, {"file_name", "tmp/rpm_fh"}}));
  mc.add(MakeCheckEnergy({{"trials_per", str(trials_per)}}));
  mc.add(MakeTune());
  mc.add(MakeEnergy({
    {"file_name", "tmp/rpm_fh_energy"},
    {"trials_per_update", "1"},
    {"trials_per_write", str(trials_per)},
    {"multistate", "true"}}));
  MonteCarlo mc2 = test_serialize(mc);
  mc2.run_until_complete();
  return mc2;
}

TEST(MonteCarlo, rpm_fh_LONG) {
  MonteCarlo mc2 = rpm_fh_test();
  FlatHistogram fh(mc2.criteria());
  const LnProbability& lnpi = fh.bias().ln_prob();
  EXPECT_NEAR(lnpi.value(0), -1.2994315780357, 0.1);
  EXPECT_NEAR(lnpi.value(1), -1.08646312498868, 0.15);
  EXPECT_NEAR(lnpi.value(2), -0.941850889679828, 0.2);
  EXPECT_NEAR(energy_av(0, mc2), 0, 1e-14);
  EXPECT_NEAR(energy_av(1, mc2), -0.939408, 0.02);
  EXPECT_NEAR(energy_av(2, mc2), -2.02625, 0.05);
}

TEST(MonteCarlo, rpm_fh_divalent_VERY_LONG) {
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
    {"alpha", str(5./15)}}));
  mc.add_to_reference(MakePotential(MakeDontVisitModel()));
//  mc.set(1, Potential(MakeModelTwoBodyFactory({MakeHardSphere(),
//                                               MakeChargeScreened()}),
//                      MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
  mc.set(MakeThermoParams({{"beta", str(1/temperature)},
     {"chemical_potential0", str(beta_mu*temperature)},
     {"chemical_potential1", str(beta_mu*temperature)}}));
  auto criteria = MakeFlatHistogram(
    MakeMacrostateNumParticles(
      Histogram({{"width", "1"}, {"max", "5"}, {"min", "0"}}),
      {{"particle_type", "0"}}),
    // MakeWangLandau({{"min_flatness", "100"}}),
    MakeTransitionMatrix({{"min_sweeps", "1000"}}));
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
  const int trials_per = 1e3;
  mc.add(MakeCriteriaUpdater({{"trials_per", str(trials_per)}}));
  mc.add(MakeCriteriaWriter({
    {"trials_per", str(trials_per)},
    {"file_name", "tmp/dival_fh_crit.txt"}}));
  mc.add(MakeLogAndMovie({{"trials_per", str(trials_per)}, {"file_name", "tmp/dival_fh"}}));
  mc.add(MakeCheckEnergy({{"trials_per", str(trials_per)}, {"tolerance", str(1e-4)}}));
  mc.add(MakeTune());
  mc.add(MakeCheckNetCharge({{"trials_per", str(trials_per)}}));
  const int en_index = mc.num_analyzers();
  mc.add(MakeEnergy({
    {"file_name", "tmp/dival_fh_energy"},
    {"trials_per_update", "1"},
    {"trials_per_write", str(trials_per)},
    {"multistate", "true"}}));
  mc.run_until_complete();
  // mc.attempt(1e7);

  const LnProbability lnpi = FlatHistogram(mc.criteria()).bias().ln_prob();
  EXPECT_NEAR(lnpi.value(0), -6.7005955776549158, 0.09);
  EXPECT_NEAR(lnpi.value(1), -3.6523345299136007, 0.07);
  EXPECT_NEAR(lnpi.value(2), -2.1178631459398805, 0.05);
  EXPECT_NEAR(lnpi.value(3), -1.3652342629553453, 0.04);
  EXPECT_NEAR(lnpi.value(4), -1.1336431696116527, 0.03);
  EXPECT_NEAR(lnpi.value(5), -1.2896341247626120, 0.03);
//  EXPECT_NEAR(lnpi.value(0), -6.6615, 0.06);
//  EXPECT_NEAR(lnpi.value(1), -3.6256, 0.06);
//  EXPECT_NEAR(lnpi.value(2), -2.1046, 0.06);
//  EXPECT_NEAR(lnpi.value(3), -1.3685, 0.06);
//  EXPECT_NEAR(lnpi.value(4), -1.1371, 0.06);
//  EXPECT_NEAR(lnpi.value(5), -1.2911, 0.06);
  const std::vector<std::shared_ptr<Analyze> >& en =
    mc.analyzers()[en_index]->analyzers();
  EXPECT_NEAR(en[0]->accumulator().average(), 0, 1e-14);
  EXPECT_NEAR(en[1]->accumulator().average(), -1.3278876302141585, 0.05);
  EXPECT_NEAR(en[2]->accumulator().average(), -3.0162868737745732, 0.05);
  EXPECT_NEAR(en[3]->accumulator().average(), -4.8648645814174927, 0.055);
  EXPECT_NEAR(en[4]->accumulator().average(), -6.8089768188067694, 0.065);
  EXPECT_NEAR(en[5]->accumulator().average(), -8.8377616317395002, 0.07);
//  EXPECT_NEAR(en[0]->accumulator().average(), 0, 1e-14);
//  EXPECT_NEAR(en[1]->accumulator().average(), -1.30701, 0.03);
//  EXPECT_NEAR(en[2]->accumulator().average(), -2.98115, 0.03);
//  EXPECT_NEAR(en[3]->accumulator().average(), -4.85254, 0.05);
//  EXPECT_NEAR(en[4]->accumulator().average(), -6.80956, 0.12);
//  EXPECT_NEAR(en[5]->accumulator().average(), -8.85025, 0.16);
}

MonteCarlo nvtw(const int num) {
  const int num_prod = 1e3;
  const int num_equil = 1e3;
  const std::string trials_per = "1e3";
  MonteCarlo mc;
  //mc.set(MakeRandomMT19937({{"seed", "1633373249"}}));
  mc.add(MakeConfiguration({{"cubic_box_length", "8"},
                            {"particle_type0", "../forcefield/lj.fstprt"}}));
  mc.add(MakePotential(MakeLennardJones()));
  mc.add(MakePotential(MakeLongRangeCorrections()));
  mc.set(MakeThermoParams({{"beta", str(1./1.5)}, {"chemical_potential", "-2.352321"}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
//  mc.add(MakeTrialGrow({{{"translate", "true"}, {"particle_type", "0"}, {"tunable_param", "1"}, {"num_steps", "4"}, {"reference_index", "0"}}}));
  mc.add(MakeTrialAdd({{"particle_type", "0"}, {"weight", "4"}}));
  mc.add(MakeTune());
  mc.add(MakeCheckEnergy({{"trials_per", str(trials_per)}}));
  mc.add(MakeLogAndMovie({{"trials_per", str(trials_per)}, {"file_name", "tmp/lj_fh"}}));
  mc.run(MakeRun({{"until_num_particles", str(num)}}));
  mc.run(MakeRemoveTrial({{"name", "TrialAdd"}}));
  mc.attempt((num+1)*num_equil);
  mc.run(MakeRemoveModify({{"name", "Tune"}}));
//  mc.add(MakeTrialGrow({{{"transfer", "true"}, {"particle_type", "0"}, {"weight", "4"}, {"num_steps", "4"}, {"reference_index", "0"}}}));
  mc.add(MakeTrialTransfer({{"particle_type", "0"}, {"weight", "4"}}));
  mc.set(MakeFlatHistogram(
    MakeMacrostateNumParticles(
      Histogram({{"width", "1"}, {"max", str(num)}, {"min", str(num)}})),
  MakeTransitionMatrix({{"min_sweeps", "1e8"}})));
  mc.add(MakeEnergy({
    {"file_name", "tmp/ljen" + str(num) + ".txt"},
    {"trials_per_update", "1"},
    {"trials_per_write", str(trials_per)}}));
  mc.add(MakeCriteriaWriter({
    {"trials_per", trials_per},
    {"file_name", "tmp/ljcrit" + str(num) + ".txt"}}));
  mc.attempt((num+1)*num_prod);
  return mc;
}

TEST(MonteCarlo, nvtw) {
  const int min = 1, max = 5;
  std::vector<std::vector<std::vector<Accumulator> > > data;
  for (int num = min; num <= max; ++num) {
    MonteCarlo mc = nvtw(num);
    FlatHistogram fh = FlatHistogram(mc.criteria());
    std::stringstream ss;
    fh.bias().serialize(ss);
    TransitionMatrix tm(ss);
    //INFO(feasst_str(tm.collection().matrix()));
    data.push_back(tm.collection().matrix());
  }
  CollectionMatrix cm(data);
  LnProbability lnp;
  lnp.resize(data.size());
  cm.compute_ln_prob(&lnp);
  //INFO(feasst_str(lnp.values()));
  EXPECT_NEAR(lnp.value(0), -14.03737335832180, 0.25);
  EXPECT_NEAR(lnp.value(1), -10.05031209165520, 0.25);
  EXPECT_NEAR(lnp.value(2), -6.458920624988570, 0.25);
  EXPECT_NEAR(lnp.value(3), -3.145637424988510, 0.15);
  EXPECT_NEAR(lnp.value(4), -0.045677458321876, 0.01);
}

}  // namespace feasst
