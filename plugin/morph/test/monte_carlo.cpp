#include "utils/test/utils.h"
#include "monte_carlo/test/monte_carlo_utils.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/configuration.h"
#include "system/include/lennard_jones.h"
#include "system/include/long_range_corrections.h"
#include "system/include/dont_visit_model.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/run.h"
#include "monte_carlo/include/remove.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/trial_add.h"
#include "flat_histogram/include/flat_histogram.h"
#include "flat_histogram/include/transition_matrix.h"
#include "flat_histogram/include/ln_probability.h"
#include "flat_histogram/include/wang_landau.h"
#include "steppers/include/criteria_writer.h"
#include "steppers/include/energy.h"
#include "steppers/include/criteria_updater.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/tune.h"
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
  mc.run(MakeRemove({{"name", "TrialAdd"}}));
  mc.add(MakeTrialAdd({{"particle_type", "1"}}));
  mc.run(MakeRun({{"until_num_particles", "4"}}));
  mc.run(MakeRemove({{"name", "TrialAdd"}}));
  EXPECT_EQ(mc.configuration().num_particles_of_type(0), 2);
  EXPECT_EQ(mc.configuration().num_particles_of_type(1), 2);
  mc.add(std::make_shared<TrialMorph>(argtype({{"particle_type0", "1"},
                         {"particle_type_morph0", "0"}, {"reference_index", "0"}})));
  mc.add(std::make_shared<TrialMorph>(argtype({{"particle_type0", "0"},
                         {"particle_type_morph0", "1"}, {"reference_index", "0"}})));
//  mc.add(MakeLogAndMovie({{"trials_per_write", str(1e2)}, {"output_file", "tmp/growth"}}));
  mc.add(MakeCheckEnergy({{"trials_per_update", str(1e2)}}));
  mc.add(MakeTune());
  mc.attempt(1e3);
//  EXPECT_EQ(mc.configuration().num_particles_of_type(0), 4);
//  EXPECT_EQ(mc.configuration().num_particles_of_type(1), 0);
}

TEST(MonteCarlo, TrialMorph) {
  System system;
  system.add(MakeConfiguration({{"cubic_side_length", "8"},
    {"particle_type0", "../particle/lj.fstprt"},
    {"particle_type1", "../particle/lj.fstprt"}}));
  system.add(MakePotential(MakeLennardJones()));
  system.add(MakePotential(MakeLongRangeCorrections()));
  system.add_to_reference(MakePotential(MakeDontVisitModel()));
  test_morph(system);
}

TEST(MonteCarlo, TrialMorphCO2N2) {
  System system;
  system.add(MakeConfiguration({{"cubic_side_length", "30"},
    {"particle_type0", "../particle/co2.fstprt"},
    {"particle_type1", "../particle/n2.fstprt"}}));
  system.add(MakePotential(MakeLennardJones()));
  system.add(MakePotential(MakeLongRangeCorrections()));
  system.add_to_reference(MakePotential(MakeDontVisitModel()));
  test_morph(system);
}

std::unique_ptr<MonteCarlo> test_morph_expanded_lj(
  const std::vector<std::vector<int> > grow_sequence,
  const int max = 5) {
  auto mc = std::make_unique<MonteCarlo>();
  mc->add(MakeConfiguration({{"cubic_side_length", "8"},
                            {"particle_type0", "../particle/lj.fstprt"},
                            {"particle_type1", "../particle/atom.fstprt"},
                            {"sigma1", "0.5"},
                            {"cutoff1", "1"},
                            {"add_particles_of_type0", "1"}}));
  mc->add(MakePotential(MakeLennardJones()));
  mc->add(MakePotential(MakeLongRangeCorrections()));
  mc->add_to_reference(MakePotential(MakeDontVisitModel()));
  const double num_parts_in_grow = static_cast<double>(grow_sequence[0].size());
  //INFO(str(num_parts_in_grow/grow_sequence.size()));
  mc->set(MakeThermoParams({
    {"beta", str(1./1.5)},
    {"chemical_potential0", "-2.352321"},
    {"chemical_potential1", "-2.352321"}}));
  auto criteria = MakeFlatHistogram(
    MakeMacrostateMorph(grow_sequence,
      Histogram({{"width", str(num_parts_in_grow/grow_sequence.size())},
                 {"max", str(max)}, {"min", "1"}})),
    MakeTransitionMatrix({{"min_sweeps", "1000"}}));
  mc->set(criteria);
  mc->add(MakeTrialTranslate({{"weight", "0.25"}, {"tunable_param", "1."}}));
  mc->add(MakeTrialMorphExpanded(grow_sequence,
    {{"reference_index", "0"}}));//, {"shift", str(-1*num_parts_in_grow)}}));
  const std::string trials_per = str(int(1e3));
//  mc->add(MakeLogAndMovie({{"trials_per_write", trials_per}, {"output_file", "tmp/grow_fh"}}));
  mc->add(MakeCheckEnergy({{"trials_per_update", trials_per}}));
  mc->add(MakeTune());
  mc->add(MakeCriteriaUpdater({{"trials_per_update", trials_per}}));
  mc->add(MakeCriteriaWriter({
    {"trials_per_write", trials_per},
    {"output_file", "tmp/grow_fh_crit.txt"}}));
  mc->add(MakeEnergy({
    {"output_file", "tmp/grow_fh_energy"},
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
  auto mc = test_morph_expanded_lj({{1, 1}, {0, 0}}, 5);
  mc->run_until_complete();
  //INFO(FlatHistogram(mc->criteria()).write());
  const LnProbability lnpi = FlatHistogram().flat_histogram(mc->criteria())->bias().ln_prob().reduce(2);
  EXPECT_NEAR(lnpi.value(0), -13.9933350923078, 0.04);
  EXPECT_NEAR(lnpi.value(1), -6.41488235897456, 0.04);
  EXPECT_NEAR(lnpi.value(2), -0.00163919230786818, 0.005);
  EXPECT_NEAR(mc->analyzers().back()->analyzers()[0]->accumulator().average(),
              -0.000605740233333333, 1e-8);
  EXPECT_NEAR(mc->analyzers().back()->analyzers()[2]->accumulator().average(),
              -0.089928316, 0.002);
  EXPECT_NEAR(mc->analyzers().back()->analyzers()[4]->accumulator().average(),
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

std::unique_ptr<MonteCarlo> test_morph_expanded(const std::string trials_per) {
  auto mc = std::make_unique<MonteCarlo>();
  mc->add(MakeConfiguration({{"cubic_side_length", "8"},
                            {"particle_type0", "../particle/lj.fstprt"},
                            {"particle_type1", "../particle/lj.fstprt"},
                            {"particle_type2", "../particle/lj.fstprt"},
                            {"particle_type3", "../particle/lj.fstprt"},
                            {"add_particles_of_type0", "1"}}));
  mc->add(MakePotential(MakeLennardJones()));
  const std::vector<std::vector<int> > grow_sequence = {{1}, {2}, {3}, {0}};
  mc->set(MakeThermoParams({
      {"beta", str(1./1.5)},
      {"chemical_potential0", "-2.352321"},
      {"chemical_potential1", "-2"},
      {"chemical_potential2", "-2.1"},
      {"chemical_potential3", "-2.2"}}));
  mc->set(MakeFlatHistogram(
    MakeMacrostateMorph(
      grow_sequence,
      Histogram({{"width", str(1./grow_sequence.size())}, {"max", "5"}, {"min", "1"}})),
    // MakeWangLandau({{"min_flatness", "25"}}),
    MakeTransitionMatrix({{"min_sweeps", "10"}})));
  mc->add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc->add(MakeTrialMorphExpanded(grow_sequence));
//  mc->add(MakeLogAndMovie({{"trials_per_write", trials_per}, {"output_file", "tmp/growth"}}));
  mc->add(MakeCheckEnergy({{"trials_per_update", trials_per}}));
  mc->add(MakeTune());
  mc->add(MakeCriteriaUpdater({{"trials_per_write", trials_per}}));
  mc->add(MakeCriteriaWriter({{"trials_per_write", trials_per}, {"output_file", "tmp/growth_crit.txt"}}));
  mc->add(MakeEnergy({{"output_file", "tmp/growth_energy"},
    {"trials_per_update", "1"},
    {"trials_per_write", str(trials_per)},
    {"multistate", "true"}}));
  return mc;
}

TEST(MonteCarlo, TrialMorphExpanded) {
  auto mc = test_morph_expanded(str(int(1e3)));
  auto mc2 = test_serialize_unique(*mc);
  EXPECT_EQ(2, mc->trials().num());
  mc2->attempt(1e3);
}

TEST(MonteCarlo, TrialMorphExpanded_LONG) {
  auto mc = test_morph_expanded(str(int(1e4)));
  mc->run_until_complete();
  // INFO(mc->criteria().write());

  std::stringstream ss;
  mc->criteria().serialize(ss);
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
  EXPECT_NEAR(energy_av2(4, *mc), -0.030574223333333334, 0.03);
  EXPECT_NEAR(energy_av2(8, *mc), -0.089928316, 0.05);
  EXPECT_NEAR(energy_av2(12, *mc), -0.1784570533333333, 0.06);
  EXPECT_NEAR(energy_av2(16, *mc), -0.29619201333333334, 0.14);
}

TEST(MonteCarlo, TrialMorphExpandedBinary_LONG) {
  MonteCarlo mc;
  mc.set(MakeRandomMT19937({{"seed", "1234"}}));
  mc.add(MakeConfiguration({{"cubic_side_length", "8"},
                            {"particle_type0", "../particle/lj.fstprt"},
                            {"particle_type1", "../particle/lj.fstprt"},
                            {"particle_type2", "../particle/lj.fstprt"},
                            {"particle_type3", "../particle/lj.fstprt"}}));
  mc.add(MakePotential(MakeLennardJones()));
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
//  mc.add(MakeLogAndMovie({{"trials_per_write", trials_per}, {"output_file", "tmp/growth"}}));
  mc.add(MakeCheckEnergy({{"trials_per_update", trials_per}}));
  mc.add(MakeTune());
  mc.add(MakeCriteriaUpdater({{"trials_per_write", trials_per}}));
  mc.add(MakeCriteriaWriter({{"trials_per_write", trials_per}, {"output_file", "tmp/growth_crit.txt"}}));
  mc.add(MakeEnergy({{"output_file", "tmp/growth_energy"},
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

TEST(MonteCarlo, morphrxn) {
  const int tpi = 1e0;
  const int trials = tpi*5e1;
  const std::string tpis = str(tpi);
  auto mc = MakeMonteCarlo({{
    {"RandomMT19937", {{"seed", "123"}}},
    {"Configuration", {{"cubic_side_length", "8"},
                       {"particle_type0", "../plugin/morph/particle/trimer_to_tetramer.fstprt"}, // h2o
                       {"particle_type1", "../particle/tetramer.fstprt"},                        // h30
                       {"particle_type2", "../plugin/morph/particle/dimer_to_monomer.fstprt"},    // dimer should go to monomer
                       {"particle_type3", "../plugin/morph/particle/monomer.fstprt"},
                       //{"particle_type0", "../particle/lj.fstprt"},
                       //{"particle_type1", "../particle/lj.fstprt"},
                       //{"particle_type2", "../particle/lj.fstprt"},
                       //{"particle_type3", "../particle/lj.fstprt"},
                       {"add_particles_of_type2", "1"}}},
    {"Potential", {{"Model", "LennardJones"}}},
    {"RefPotential", {{"VisitModel", "DontVisitModel"}}},
    {"ThermoParams", {{"beta", "1"}, {"chemical_potential0", "1."}, {"chemical_potential1", "1."},
                                     {"chemical_potential2", "1."}, {"chemical_potential3", "1."}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"weight_per_number_fraction", str(1./8.)}, {"particle_type", "0"}}},
    {"TrialTranslate", {{"weight_per_number_fraction", str(1./8.)}, {"particle_type", "1"}}},
    {"TrialTranslate", {{"weight_per_number_fraction", str(1./8.)}, {"particle_type", "2"}}},
    {"TrialTranslate", {{"weight_per_number_fraction", str(1./8.)}, {"particle_type", "3"}}},
    //{"TrialTranslate", {{"weight_per_number_fraction", "1."}, {"particle_type", "1"}, {"number_fraction_exclude_type0", "2"}, {"number_fraction_exclude_type1", "3"}}},
    {"TrialParticlePivot", {{"weight_per_number_fraction", str(1./8.)}, {"particle_type", "0"}}},
    {"TrialParticlePivot", {{"weight_per_number_fraction", str(1./8.)}, {"particle_type", "1"}}},
    {"TrialParticlePivot", {{"weight_per_number_fraction", str(1./8.)}, {"particle_type", "2"}}},
    {"TrialParticlePivot", {{"weight_per_number_fraction", str(1./8.)}, {"particle_type", "3"}}},
    {"TrialAdd", {{"particle_type", "0"}}},
    {"Run", {{"until_num_particles", "100"}, {"particle_type", "0"}}},
    {"Remove", {{"name", "TrialAdd"}}},
    {"TrialMorph", {{"weight", "0.1"}, {"reference_index", "0"},
                    {"particle_type0", "0"}, {"particle_type_morph0", "1"},
                    {"particle_type1", "2"}, {"particle_type_morph1", "3"}}},
    {"TrialMorph", {{"weight", "0.1"}, {"reference_index", "0"},
                    {"particle_type0", "1"}, {"particle_type_morph0", "0"},
                    {"particle_type1", "3"}, {"particle_type_morph1", "2"}}},
    {"TrialMorph", {{"weight", "0.1"}, {"reference_index", "0"},
                    {"particle_type0", "0"}, {"particle_type_morph0", "1"},
                    {"particle_type1", "1"}, {"particle_type_morph1", "0"}}},
    {"CheckEnergy", {{"trials_per_update", tpis}, {"decimal_places", "8"}}},
    //{"Log", {{"trials_per_write", tpis}, {"output_file", "tmp/lj.csv"}}},
    //{"Movie", {{"trials_per_write", tpis}, {"output_file", "tmp/lj.xyz"}}},
    //{"NumParticles", {{"trials_per_write", tpis}, {"output_file", "tmp/ljn.csv"}, {"particle_type", "2"}}},
    //{"ProfileCPU", {{"trials_per_write", tpis}, {"output_file", "tmp/ljp.csv"}, {"trials_per_update", "1e3"}}},
    {"Run", {{"num_trials", str(trials)}}},
  }}, true);
}

TEST(MonteCarlo, octane_01_fh_VERY_LONG) {
  const std::string tpc = "1e4";
  auto mc = MakeMonteCarlo({{
    {"Configuration", {{"cubic_side_length", "45"}, {"particle_type0", "../particle/n-octane.fstprt"}}},
    {"Potential", {{"Model", "LennardJones"}}},
    {"Potential", {{"Model", "LennardJones"}, {"VisitModel", "VisitModelIntraMap"}, {"exclude_bonds", "true"}, {"exclude_angles", "true"}, {"exclude_dihedrals", "true"}}},
    {"RefPotential", {{"VisitModel", "DontVisitModel"}}},
    {"ThermoParams", {{"beta", "0.3436353001220744"}, {"chemical_potential0", "-17.460371498121805"}}},
    {"FlatHistogram", {{"Macrostate", "MacrostateNumParticles"}, {"width", "1"}, {"max", "1"}, {"min", "0"},
      {"Bias", "TransitionMatrix"}, {"min_sweeps", "1e2"}}},
    {"TrialGrowFile", {{"grow_file", "../plugin/chain/test/data/trappe_grow_grand_canonical.txt"}}},
    {"Energy", {{"trials_per_write", tpc}, {"output_file", "tmp/oct_en.txt"}, {"multistate", "true"}}},
    {"Log", {{"trials_per_write", tpc}, {"output_file", "tmp/oct.txt"}}},
    {"Movie", {{"trials_per_write", tpc}, {"output_file", "tmp/oct.xyz"}}},
    {"CheckEnergy", {{"trials_per_update", tpc}, {"tolerance", str(1e-9)}}},
    {"CriteriaUpdater", {{"trials_per_update", tpc}}},
    {"CriteriaWriter", {{"trials_per_write", tpc}, {"output_file", "tmp/oct.csv"}}},
    {"Tune", {{}}},
    {"Run", {{"until", "complete"}}},
  }}, true);
  const LnProbability lnpi = FlatHistogram().flat_histogram(mc->criteria())->bias().ln_prob();
  EXPECT_NEAR(lnpi.delta(1), 5.8699, 5*0.00749341);
  const std::vector<std::shared_ptr<Analyze> >& en = mc->analyzers()[0]->analyzers();
  EXPECT_NEAR(en[0]->accumulator().average(), 0, 1e-13);
  EXPECT_NEAR(en[1]->accumulator().average(), 20.955607316224864, 5*0.061054957999540201);
  EXPECT_GT(mc->trial(0).num_success(), 5e4);
}

}  // namespace feasst
