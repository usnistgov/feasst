#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "system/include/lennard_jones.h"
#include "system/include/long_range_corrections.h"
#include "system/include/dont_visit_model.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/run.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/trial_add.h"
#include "flat_histogram/include/flat_histogram.h"
#include "flat_histogram/include/transition_matrix.h"
#include "flat_histogram/include/wang_landau.h"
#include "steppers/include/criteria_writer.h"
#include "steppers/include/energy.h"
#include "steppers/include/criteria_updater.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/tune.h"
#include "steppers/include/log_and_movie.h"
#include "morph/include/trial_morph.h"
#include "morph/include/trial_morph_expanded.h"
#include "morph/include/macrostate_morph.h"
#include "charge/include/utils.h"
#include "charge/test/charge_utils.h"

namespace feasst {

double energy_av2(const int macro, const MonteCarlo& mc) {
  return mc.analyzers().back()->analyzers()[macro]->accumulator().average();
}

void test_morph(const System& system) {
  MonteCarlo mc;
  mc.set(system);
  mc.set(MakeThermoParams({{"beta", "1"}, {"chemical_potential0", "1."},
                                        {"chemical_potential1", "1."}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc.add(MakeTrialAdd({{"particle_type", "0"}}));
  mc.run(MakeRun({{"until_num_particles", "2"}}));
  mc.run(MakeRemoveTrial({{"name", "TrialAdd"}}));
  mc.add(MakeTrialAdd({{"particle_type", "1"}}));
  mc.run(MakeRun({{"until_num_particles", "4"}}));
  mc.run(MakeRemoveTrial({{"name", "TrialAdd"}}));
  EXPECT_EQ(mc.configuration().num_particles_of_type(0), 2);
  EXPECT_EQ(mc.configuration().num_particles_of_type(1), 2);
  mc.add(MakeTrialMorph({{"particle_type0", "1"},
                                 {"particle_type_morph0", "0"}}));
  mc.add(MakeLogAndMovie({{"trials_per", str(1e2)}, {"file_name", "tmp/growth"}}));
  mc.add(MakeCheckEnergy({{"trials_per", str(1e2)}}));
  mc.add(MakeTune());
  mc.attempt(1e3);
  EXPECT_EQ(mc.configuration().num_particles_of_type(0), 4);
  EXPECT_EQ(mc.configuration().num_particles_of_type(1), 0);
}

TEST(MonteCarlo, TrialMorph) {
  System system;
  system.add(*MakeConfiguration({{"cubic_box_length", "8"}, {"particle_type0", "../forcefield/lj.fstprt"}}));
  system.add(MakePotential(MakeLennardJones()));
  system.add(MakePotential(MakeLongRangeCorrections()));
  Configuration * config = system.get_configuration();
  config->add_particle_type(install_dir() + "/forcefield/lj.fstprt", "0.5");
  config->set_model_param("sigma", 1, 0.5);
  config->set_model_param("cutoff", 1, 1.0);
  test_morph(system);
}

MonteCarlo test_morph_expanded_lj(
  const std::vector<std::vector<int> > grow_sequence,
  const int max = 5) {
  MonteCarlo mc;
  {
    auto config = MakeConfiguration({{"cubic_box_length", "8"}, {"particle_type0", "../forcefield/lj.fstprt"}});
    config->add_particle_type(install_dir() + "/forcefield/lj.fstprt", "0.5");
    config->set_model_param("sigma", 1, 0.5);
    config->set_model_param("cutoff", 1, 1.0);
    config->add_particle_of_type(0);
    mc.add(config);
    mc.add(MakePotential(MakeLennardJones()));
    mc.add(MakePotential(MakeLongRangeCorrections()));
  }
  mc.add_to_reference(MakePotential(MakeDontVisitModel()));
  const double num_parts_in_grow = static_cast<double>(grow_sequence[0].size());
  //INFO(str(num_parts_in_grow/grow_sequence.size()));
  mc.set(MakeThermoParams({
    {"beta", str(1./1.5)},
    {"chemical_potential0", "-2.352321"},
    {"chemical_potential1", "-2.352321"}}));
  auto criteria = MakeFlatHistogram(
    MakeMacrostateMorph(grow_sequence,
      Histogram({{"width", str(num_parts_in_grow/grow_sequence.size())},
                 {"max", str(max)}, {"min", "1"}})),
    MakeTransitionMatrix({{"min_sweeps", "1000"}}));
  mc.set(criteria);
  mc.add(MakeTrialTranslate({{"weight", "0.25"}, {"tunable_param", "1."}}));
  mc.add(MakeTrialMorphExpanded(grow_sequence,
    {{"reference_index", "0"}}));//, {"shift", str(-1*num_parts_in_grow)}}));
  const std::string trials_per = str(int(1e3));
  mc.add(MakeLogAndMovie({{"trials_per", trials_per}, {"file_name", "tmp/grow_fh"}}));
  mc.add(MakeCheckEnergy({{"trials_per", trials_per}}));
  mc.add(MakeTune());
  mc.add(MakeCriteriaUpdater({{"trials_per", trials_per}}));
  mc.add(MakeCriteriaWriter({
    {"trials_per", trials_per},
    {"file_name", "tmp/grow_fh_crit.txt"}}));
  mc.add(MakeEnergy({
    {"file_name", "tmp/grow_fh_energy"},
    {"trials_per_update", "1"},
    {"trials_per_write", trials_per},
    {"multistate", "true"}}));
  return mc;
}

//TEST(MonteCarlo, TrialMorphExpanded_args) {
//  try {
//    TrialMorphExpanded({{1, -1}, {0, 0}});
//    CATCH_PHRASE("cant have add and morph in same stage");
//  }
//}

TEST(MonteCarlo, TrialMorphExpanded_2_lj_LONG) {
  MonteCarlo mc = test_morph_expanded_lj({{1, 1}, {0, 0}}, 5);
  mc.run_until_complete();
  INFO(FlatHistogram(mc.criteria()).write());
  const LnProbability lnpi = FlatHistogram(mc.criteria()).bias().ln_prob().reduce(2);
  EXPECT_NEAR(lnpi.value(0), -13.9933350923078, 0.04);
  EXPECT_NEAR(lnpi.value(1), -6.41488235897456, 0.04);
  EXPECT_NEAR(lnpi.value(2), -0.00163919230786818, 0.005);
  EXPECT_NEAR(mc.analyzers().back()->analyzers()[0]->accumulator().average(),
              -0.000605740233333333, 1e-8);
  EXPECT_NEAR(mc.analyzers().back()->analyzers()[2]->accumulator().average(),
              -0.089928316, 0.002);
  EXPECT_NEAR(mc.analyzers().back()->analyzers()[4]->accumulator().average(),
              -0.29619201333333334, 0.02);
}

//TEST(MonteCarlo, TrialMorphExpanded_2_skip_lj_LONG) {
//  MonteCarlo mc = test_morph_expanded_lj({{1, -1},
//                                          {0, -1},
//                                          {-1, 1},
//                                          {-1, 0}});
//  mc.run_until_complete();
//  INFO(FlatHistogram(mc.criteria()).write());
//  const LnProbability lnpi = FlatHistogram(mc.criteria()).bias().ln_prob().reduce(4);
//  EXPECT_NEAR(lnpi.value(0), -13.9933350923078, 0.02);
//  EXPECT_NEAR(lnpi.value(1), -6.41488235897456, 0.02);
//  EXPECT_NEAR(lnpi.value(2), -0.00163919230786818, 0.005);
//  EXPECT_NEAR(mc.analyzers().back()->analyzers()[0]->accumulator().average(),
//              -0.000605740233333333, 1e-8);
//  EXPECT_NEAR(mc.analyzers().back()->analyzers()[4]->accumulator().average(),
//              -0.089928316, 0.002);
//  EXPECT_NEAR(mc.analyzers().back()->analyzers()[8]->accumulator().average(),
//              -0.29619201333333334, 0.02);
//}

//HWH need to fix macrostate
//TEST(MonteCarlo, TrialMorphExpanded_3_lj_LONG) {
//  MonteCarlo mc = test_morph_expanded_lj({{1, 1, 1}, {0, 0, 0}}, 4);
//  mc.run_until_complete();
//  INFO(FlatHistogram(mc.criteria()).write());
//  const LnProbability lnpi = FlatHistogram(mc.criteria()).bias().ln_prob().reduce(3);
//  EXPECT_NEAR(lnpi.value(0), -10.8917545445662, 0.02);
//  EXPECT_NEAR(lnpi.value(1), -0.000018611232922473, 0.02);
//  EXPECT_NEAR(mc.analyzers().back()->analyzers()[0]->accumulator().average(),
//              -0.000605740233333333, 1e-8);
//  EXPECT_NEAR(mc.analyzers().back()->analyzers()[3]->accumulator().average(),
//              -0.1784570533333333, 0.02);
//}

TEST(MonteCarlo, TrialMorph_RPM) {
  test_morph(rpm({{"alpha", str(5.6/12)}, {"kmax_squared", "38"}}));
}

MonteCarlo test_morph_expanded(const std::string trials_per) {
  MonteCarlo mc;
  //mc.set(MakeRandomMT19937({{"seed", "1234"}}));
  {
    auto config = MakeConfiguration({{"cubic_box_length", "8"}, {"particle_type0", "../forcefield/lj.fstprt"}});
    config->add_particle_type(install_dir() + "/forcefield/lj.fstprt", "0.25");
    config->set_model_param("sigma", 1, 0.25);
    config->set_model_param("cutoff", 1, 1.0);
    config->add_particle_type(install_dir() + "/forcefield/lj.fstprt", "0.5");
    config->set_model_param("sigma", 1, 0.5);
    config->set_model_param("cutoff", 1, 1.0);
    config->add_particle_type(install_dir() + "/forcefield/lj.fstprt", "0.75");
    config->set_model_param("sigma", 1, 0.75);
    config->set_model_param("cutoff", 1, 1.0);
    config->add_particle_of_type(0);
    mc.add(config);
    mc.add(MakePotential(MakeLennardJones()));
  }
  const std::vector<std::vector<int> > grow_sequence = {{1}, {2}, {3}, {0}};
  mc.set(MakeThermoParams({
      {"beta", str(1./1.5)},
      {"chemical_potential0", "-2.352321"},
      {"chemical_potential1", "-2"},
      {"chemical_potential2", "-2.1"},
      {"chemical_potential3", "-2.2"}}));
  mc.set(MakeFlatHistogram(
    MakeMacrostateMorph(
      grow_sequence,
      Histogram({{"width", str(1./grow_sequence.size())}, {"max", "5"}, {"min", "1"}})),
    // MakeWangLandau({{"min_flatness", "25"}}),
    MakeTransitionMatrix({{"min_sweeps", "10"}})));
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc.add(MakeTrialMorphExpanded(grow_sequence));
  mc.add(MakeLogAndMovie({{"trials_per", trials_per}, {"file_name", "tmp/growth"}}));
  mc.add(MakeCheckEnergy({{"trials_per", trials_per}}));
  mc.add(MakeTune());
  mc.add(MakeCriteriaUpdater({{"trials_per_write", trials_per}}));
  mc.add(MakeCriteriaWriter({{"trials_per_write", trials_per}, {"file_name", "tmp/growth_crit.txt"}}));
  mc.add(MakeEnergy({{"file_name", "tmp/growth_energy"},
    {"trials_per_update", "1"},
    {"trials_per_write", str(trials_per)},
    {"multistate", "true"}}));
  return mc;
}

TEST(MonteCarlo, TrialMorphExpanded) {
  MonteCarlo mc = test_morph_expanded(str(int(1e3)));
  MonteCarlo mc2 = test_serialize(mc);
  EXPECT_EQ(2, mc.trials().num());
  mc2.attempt(1e3);
}

TEST(MonteCarlo, TrialMorphExpanded_LONG) {
  MonteCarlo mc = test_morph_expanded(str(int(1e4)));
  mc.run_until_complete();
  // INFO(mc.criteria().write());

  std::stringstream ss;
  mc.criteria().serialize(ss);
  FlatHistogram fh(ss);
  LnProbability lnpi = fh.bias().ln_prob().reduce(4);
  lnpi.normalize();
//  INFO(feasst_str(lnpi.values()));

  // Note: LRCs are disabled but these values assume LRCs
  // copied from flat_histogram/test/monte_carlo.cpp
  EXPECT_NEAR(lnpi.value(0), -14.037373358321800000, 0.75);
  EXPECT_NEAR(lnpi.value(1), -10.050312091655200000, 0.6);
  EXPECT_NEAR(lnpi.value(2), -6.458920624988570000, 0.55);
  EXPECT_NEAR(lnpi.value(3), -3.145637424988510000, 0.55);
  EXPECT_NEAR(lnpi.value(4), -0.045677458321876000, 0.55);

//  EXPECT_NEAR(energy_av2(0, mc), -0.000605740233333333, 1e-8);
  EXPECT_NEAR(energy_av2(4, mc), -0.030574223333333334, 0.03);
  EXPECT_NEAR(energy_av2(8, mc), -0.089928316, 0.05);
  EXPECT_NEAR(energy_av2(12, mc), -0.1784570533333333, 0.06);
  EXPECT_NEAR(energy_av2(16, mc), -0.29619201333333334, 0.14);
}

TEST(MonteCarlo, TrialMorphExpandedBinary_LONG) {
  MonteCarlo mc;
  mc.set(MakeRandomMT19937({{"seed", "1234"}}));
  {
    auto config = MakeConfiguration({{"cubic_box_length", "8"}, {"particle_type0", "../forcefield/lj.fstprt"}});
    config->add_particle_type(install_dir() + "/forcefield/lj.fstprt", "b");
    config->add_particle_type(install_dir() + "/forcefield/lj.fstprt", "0.5");
    config->set_model_param("sigma", 1, 0.5);
    config->set_model_param("cutoff", 1, 1.0);
    config->add_particle_type(install_dir() + "/forcefield/lj.fstprt", "0.5b");
    config->set_model_param("sigma", 1, 0.5);
    config->set_model_param("cutoff", 1, 1.0);
    mc.add(config);
    mc.add(MakePotential(MakeLennardJones()));
  }
  mc.add_to_reference(MakePotential(MakeDontVisitModel()));
  const std::vector<std::vector<int> > grow_sequence = {{2, 3, 3}, {0, 1, 1}};
  mc.set(MakeThermoParams({
    {"beta", str(1./1.5)},
    {"chemical_potential0", "-2.352321"},
    {"chemical_potential1", "-2"},
    {"chemical_potential2", "-2.1"},
    {"chemical_potential3", "-2.2"}}));
  mc.set(MakeFlatHistogram(
    MakeMacrostateMorph(
      grow_sequence,
      Histogram({{"width", str(1./grow_sequence.size())}, {"max", "5"}, {"min", "0"}})),
    // MakeWangLandau({{"min_flatness", "25"}}),
    MakeTransitionMatrix({{"min_sweeps", "10"}})));
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc.add(MakeTrialMorphExpanded(grow_sequence, {{"reference_index", "0"}}));
  const std::string trials_per = str(int(1e3));
  mc.add(MakeLogAndMovie({{"trials_per", trials_per}, {"file_name", "tmp/growth"}}));
  mc.add(MakeCheckEnergy({{"trials_per", trials_per}}));
  mc.add(MakeTune());
  mc.add(MakeCriteriaUpdater({{"trials_per_write", trials_per}}));
  mc.add(MakeCriteriaWriter({{"trials_per_write", trials_per}, {"file_name", "tmp/growth_crit.txt"}}));
  mc.add(MakeEnergy({{"file_name", "tmp/growth_energy"},
    {"trials_per_update", "1"},
    {"trials_per_write", str(trials_per)},
    {"multistate", "true"}}));
  for (int i = 0; i < 1e6; ++i) {
    mc.attempt(1);
    ASSERT(2*mc.configuration().num_particles_of_type(0) ==
             mc.configuration().num_particles_of_type(1), "er");
    ASSERT(2*mc.configuration().num_particles_of_type(2) ==
             mc.configuration().num_particles_of_type(3), "er");
  }
}

}  // namespace feasst
