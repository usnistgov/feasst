#include <iostream>
#include <memory>
#include "utils/test/utils.h"
#include "utils/include/io.h"
#include "math/include/accumulator.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/domain.h"
#include "system/include/visit_model_inner.h"
#include "system/include/thermo_params.h"
#include "system/include/potential.h"
#include "system/include/ideal_gas.h"
#include "system/include/hard_sphere.h"
#include "system/include/lennard_jones.h"
#include "system/include/long_range_corrections.h"
#include "system/include/visit_model_intra_map.h"
#include "system/include/visit_model_cutoff_outer.h"
#include "system/include/dont_visit_model.h"
#include "models/include/square_well.h"
#include "models/include/lennard_jones_force_shift.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/test/monte_carlo_utils.h"
#include "monte_carlo/include/trial_rotate.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/trial_add.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/run.h"
#include "monte_carlo/include/remove.h"
#include "monte_carlo/include/trial_compute_move.h"
#include "monte_carlo/include/trial_select_dihedral.h"
#include "monte_carlo/include/perturb_dihedral.h"
#include "monte_carlo/include/rosenbluth.h"
#include "steppers/include/log.h"
#include "steppers/include/tune.h"
#include "steppers/include/movie.h"
#include "steppers/include/wrap_particles.h"
#include "steppers/include/tune.h"
#include "steppers/include/check.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/check_physicality.h"
#include "steppers/include/check_properties.h"
#include "steppers/include/chirality_2d.h"
#include "steppers/include/energy.h"
#include "steppers/include/seek_analyze.h"
#include "cluster/include/energy_map_all.h"
#include "cluster/include/energy_map_neighbor.h"
#include "cluster/include/energy_map_neighbor_criteria.h"
#include "cluster/include/trial_avb2.h"
#include "chain/include/trial_pivot.h"
#include "chain/include/trial_particle_pivot.h"
#include "chain/include/trial_crankshaft.h"
#include "chain/include/trial_crankshaft_small.h"
#include "chain/include/trial_grow.h"
#include "chain/include/trial_grow_linear.h"
#include "chain/include/analyze_bonds.h"
#include "chain/test/system_chain.h"
#include "mayer/include/mayer_sampling.h"

namespace feasst {

TEST(MonteCarlo, chain) {
  const int trials_per = 1e0;
  MonteCarlo mc;
  mc.set(MakeRandomMT19937({{"seed", "1610724059"}}));
  mc.set(chain_system());
  mc.set(MakeThermoParams({{"beta", "1"}, {"chemical_potential", "1."}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialAdd({{"particle_type", "0"}}));
  mc.run(MakeRun({{"until_num_particles", "1"}}));
  mc.run(MakeRemove({{"name", "TrialAdd"}}));
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc.add(MakeTrialRotate({{"weight", "1."}, {"tunable_param", "20."}}));
  mc.add(MakeTrialPivot({
    {"weight", "1."},
    {"tunable_param", "20."},
    {"max_length", "30"},
    {"reference_index", "0"},
    {"num_steps", "2"}}));
//  // HWH Reptate finalize needs to be made to work with cell list.
//  mc.add(MakeTrialReptate({
//    {"weight", "1."},
//    {"max_length", "1"},
////    {"reference_index", "0"},
////    {"num_steps", "2"},
//  }));
  mc.add(MakeTrialCrankshaft({
    {"weight", "1."},
    {"tunable_param", "25."},
    {"max_length", "5."},
    {"reference_index", "0"},
    {"num_steps", "2"}}));
  mc.add(MakeTrialGrowLinear(
    MakeTrialComputeMove(),
    {
//      {"weight", "0.1"},
      {"particle_type", "0"},
      {"num_steps", "3"},
      {"reference_index", "0"},
    }));
  mc.add(MakeLog({
    {"trials_per_write", str(trials_per)},
    {"output_file", "tmp/chainlog.txt"}}));
  mc.add(MakeMovie({
    {"trials_per_write", str(trials_per)},
    {"output_file", "tmp/chain10movie.xyz"}}));
  mc.add(MakeCheckEnergy({
    {"trials_per_update", str(trials_per)},
    {"tolerance", "1e-10"}}));
  mc.add(MakeTune());
  mc.attempt(3e2);

  auto mc2 = test_serialize_unique(mc);
  EXPECT_EQ(mc2->analyzers().size(), 2);
  mc2->attempt(3e2);
}

// HWH this test is known to fail infrequently
// HWH sometimes it happens with valgrind and serialization, which is why the test is labelled long
TEST(MonteCarlo, TrialGrow_LONG) {
  for (const std::string particle : {"lj", "spce"}) {
    INFO(particle);
    double box_length = 8.;
    std::string data = "../particle/dimer.txt";
    if (particle == "spce") {
      box_length=20;
      data = "../particle/spce.txt";
    }
    MonteCarlo mc;
    //mc.set(MakeRandomMT19937({{"seed", "1635356012"}}));
    mc.add(MakeConfiguration({{"cubic_side_length", str(box_length)},
                              {"particle_type0", data}}));
    mc.add(MakePotential(MakeLennardJones()));
    mc.add(MakePotential(MakeLongRangeCorrections()));
    mc.add_to_reference(MakePotential(MakeLennardJones()));
    mc.set(MakeThermoParams({{"beta", "1.2"}, {"chemical_potential", "1."}}));
    mc.set(MakeMetropolis());
    mc.add(MakeTrialAdd({{"particle_type", "0"}}));
    mc.run(MakeRun({{"until_num_particles", "3"}}));
    mc.run(MakeRemove({{"name", "TrialAdd"}}));
    mc.set(MakeThermoParams({{"beta", "1.2"}, {"chemical_potential", "-700"}}));
    std::vector<argtype> grow_args;
    if (particle == "lj") {
      grow_args = {
      {{"default_num_steps", "4"}, {"default_reference_index", "0"},
       {"transfer", "true"},
       {"regrow", "true"},
       {"particle_type", "0"},
       {"site", "D1"},
       {"weight", "100"}},
      {{"bond", "true"},
       {"mobile_site", "D2"},
       {"anchor_site", "D1"},
       {"num_steps", "5"}}};
    } else if (particle == "spce") {
      grow_args = {
      {{"default_num_steps", "4"}, {"default_reference_index", "0"},
       {"transfer", "true"},
       {"regrow", "true"},
       {"particle_type", "0"},
       {"site", "0"},
       {"weight", "100"}},
      {{"bond", "true"},
       {"mobile_site", "1"},
       {"anchor_site", "0"},
       {"num_steps", "5"}},
      {{"angle", "true"},
       {"mobile_site", "2"},
       {"anchor_site", "0"},
       {"anchor_site2", "1"}}};
    }
    mc.add(MakeTrialGrow(grow_args));
    EXPECT_EQ(4, mc.trial(0).stage(0).rosenbluth().num());
    EXPECT_EQ(5, mc.trial(0).stage(1).rosenbluth().num());
    EXPECT_EQ(3, mc.trials().num());
    EXPECT_EQ(50, mc.trial(0).weight());
    EXPECT_EQ(50, mc.trial(1).weight());
    EXPECT_EQ(100, mc.trial(2).weight());
//    mc.add(MakeLogAndMovie({{"trials_per_write", str(1e0)}, {"output_file", "tmp/lj"}}));
    mc.add(MakeCheckEnergy({{"trials_per_update", str(1e0)}, {"tolerance", str(1e-9)}}));
    mc.add(MakeTune());
    EXPECT_EQ(3, mc.trials().num());
    EXPECT_TRUE(mc.trial(0).stage(0).trial_select().is_ghost());   // add
    EXPECT_FALSE(mc.trial(1).stage(0).trial_select().is_ghost());  // remove
    EXPECT_FALSE(mc.trial(2).stage(0).trial_select().is_ghost());  // regrow
    mc.add(MakeMovie(
     {{"trials_per_write", "1"},
      {"output_file", "tmp/grow.xyz"}}));
    for (int trial = 0; trial < 2e1; ++trial) {
      mc.attempt(1);
      //mc.configuration().check();
    }
    EXPECT_LT(mc.configuration().num_particles(), 4);
    mc.set(MakeThermoParams({{"beta", "1.2"}, {"chemical_potential", "100"}}));
    auto mc2 = test_serialize_unique(mc);
    mc2->attempt(2e1);
    EXPECT_GE(mc2->configuration().num_particles(), 1);
    mc2->configuration().check();
    // INFO(mc.trial(1)->accept().perturbed().str());
  }
}

std::unique_ptr<MonteCarlo> cg7mab2(const std::string& data, const int num, const int trials_per = 1) {
  INFO("data " << data);
  std::unique_ptr<MonteCarlo> mc = std::make_unique<MonteCarlo>();
  //mc.set(MakeRandomMT19937({{"seed", "1633624429"}}));
  mc->add(MakeConfiguration({{"cubic_side_length", "30"},
    {"particle_type", "../plugin/chain/particle/" + data},
    {"set_cutoff_min_to_sigma", "true"}}));
  mc->add(MakePotential(MakeHardSphere()));
  if (is_found_in(data, "fullangflex")) {
    mc->add(MakePotential(MakeHardSphere(),
                     MakeVisitModelIntraMap({{"exclude_bonds", "true"}})));
  }
  mc->set(MakeThermoParams({{"beta", "1."}, {"chemical_potential", "1"}}));
  mc->set(MakeMetropolis());
  mc->add(MakeTrialAdd({{"particle_type", "0"}}));
  mc->run(MakeRun({{"until_num_particles", str(num)}}));
  mc->run(MakeRemove({{"name", "TrialAdd"}}));
  if (is_found_in(data, "fullangflex")) {
    mc->add(MakeTrialGrow({
      //{{"particle_type", "0"}, {"site", "0"}, {"weight", "4"}, {"regrow", "1"}},
      //{{"bond", "1"}, {"mobile_site", "2"}, {"anchor_site", "0"}},
      {{"particle_type", "0"}, {"weight", "4"}, {"bond", "1"}, {"mobile_site", "1"}, {"anchor_site", "0"}},
      {{"angle", "1"}, {"mobile_site", "2"}, {"anchor_site", "1"}, {"anchor_site2", "0"}},
      {{"branch", "1"}, {"mobile_site", "3"}, {"mobile_site2", "5"}, {"anchor_site", "0"}, {"anchor_site2", "1"}},
      {{"angle", "1"}, {"mobile_site", "4"}, {"anchor_site", "3"}, {"anchor_site2", "0"}},
      {{"angle", "1"}, {"mobile_site", "6"}, {"anchor_site", "5"}, {"anchor_site2", "0"}}}));
  } else {
    mc->add(MakeTrialGrow({
      //{{"particle_type", "0"}, {"site", "0"}, {"weight", "4"}, {"regrow", "1"}},
      //{{"bond", "1"}, {"mobile_site", "2"}, {"anchor_site", "0"}},
      {{"particle_type", "0"}, {"weight", "4"}, {"bond", "1"}, {"mobile_site", "2"}, {"anchor_site", "0"}},
      {{"angle", "1"}, {"mobile_site", "1"}, {"anchor_site", "2"}, {"anchor_site2", "0"}},
      {{"branch", "1"}, {"mobile_site", "4"}, {"mobile_site2", "6"}, {"anchor_site", "0"}, {"anchor_site2", "2"}},
      {{"angle", "1"}, {"mobile_site", "3"}, {"anchor_site", "4"}, {"anchor_site2", "0"}},
      {{"angle", "1"}, {"mobile_site", "5"}, {"anchor_site", "6"}, {"anchor_site2", "0"}}}));
  }
//  mc->add(MakeLogAndMovie({{"trials_per_write", str(trials_per)}, {"output_file", "tmp/" + data}}));
  mc->add(MakeCheckEnergy({{"trials_per_update", str(trials_per)}, {"tolerance", str(1e-9)}}));
  mc->add(MakeTune());
  return mc;
}

TEST(MonteCarlo, cg7mab2) {
  for (const std::string data : {
      "cg7mab2.txt",
      "cg7mab2closed.txt",
      "cg7mab2extended.txt",
      "cg7mab2flex.txt",
      "cg7mab2fullangflex.txt",
    }) {
    cg7mab2(data, 1)->attempt(1e2);
  }
}

TEST(MonteCarlo, cg7mab2_LONG) {
  for (const std::string data : {
      "cg7mab2.txt",
      "cg7mab2closed.txt",
      "cg7mab2extended.txt",
      "cg7mab2flex.txt",
      "cg7mab2fullangflex.txt",
    }) {
    cg7mab2(data, 10, 1e4)->attempt(1e6);
  }
}

TEST(System, Angles2D) {
  MonteCarlo mc;
  mc.add(MakeConfiguration({{"side_length", "6,6"},
    {"particle_type", "../plugin/chain/particle/heterotrimer2d.txt"},
    {"add_particles_of_type0", "1"}}));
  mc.add(MakePotential(MakeLennardJones()));
  mc.set(MakeThermoParams({{"beta", "1"}}));
  mc.set(MakeMetropolis());
  DEBUG(mc.criteria().current_energy());
}

std::unique_ptr<MonteCarlo> test_avb(const bool avb2, const bool avb4 = true) {
  MonteCarlo mc;
  mc.set(MakeRandomMT19937({{"seed", "time"}}));
  //mc.set(MakeRandomMT19937({{"seed", "123"}}));
  mc.add(MakeConfiguration({{"side_length", "6,6"},
    {"particle_type", "../plugin/chain/particle/heterotrimer2d.txt"}}));
  EXPECT_EQ(2, mc.configuration().dimension());
  mc.add(MakePotential(MakeLennardJones(),
    MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
    //MakeVisitModel(MakeVisitModelInner(MakeEnergyMapNeighbor()))));
  mc.set(MakeThermoParams({{"beta", "1"}, {"chemical_potential", "1"}}));
  mc.set(MakeMetropolis());
  if (!avb2 && !avb4) {
    mc.add(MakeTrialTranslate());
    mc.add(MakeTrialRotate());
  }
//  mc.add(MakeTrialRigidCluster(
//    neighbor_criteria,
//    { {"rotate_param", "50"},
//      {"translate_param", "1"}}));
//  mc.add(MakeTrialGrow({
//    {{"regrow", "true"}, {"particle_type", "0"}, {"site", "0"}},
//    {{"bond", "true"}, {"mobile_site", "1"}, {"anchor_site", "0"}},
//    {{"angle", "true"}, {"mobile_site", "2"}, {"anchor_site", "0"}, {"anchor_site2", "1"}}},
//    {{"num_steps", "1"}}));
//  mc.add(MakeTrialGrow({
//    {{"bond", "true"}, {"particle_type", "0"}, {"mobile_site", "1"}, {"anchor_site", "0"}},
//    {{"angle", "true"}, {"mobile_site", "2"}, {"anchor_site", "0"}, {"anchor_site2", "1"}}},
//    {{"num_steps", "1"}}));
//  mc.add(MakeTrialGrow({
//    {{"bond", "true"}, {"particle_type", "0"}, {"mobile_site", "2"}, {"anchor_site", "0"}},
//    {{"angle", "true"}, {"mobile_site", "1"}, {"anchor_site", "0"}, {"anchor_site2", "2"}}},
//    {{"num_steps", "1"}}));
  mc.add(MakeNeighborCriteria({{"maximum_distance", "0.75"}, {"minimum_distance", "0.5"},
    {"site_type0", "1"}, {"site_type1", "1"}}));
  if (avb2) mc.add(MakeTrialGrow({
    {{"default_num_steps", "1"},
     {"regrow_avb2", "true"}, {"particle_type", "0"}, {"site", "1"}, {"neighbor_index", "0"}, {"target_particle_type", "0"}, {"target_site", "1"}},
    {{"bond", "true"}, {"mobile_site", "0"}, {"anchor_site", "1"}},
    {{"angle", "true"}, {"mobile_site", "2"}, {"anchor_site", "0"}, {"anchor_site2", "1"}}}));
  if (avb4) mc.add(MakeTrialGrow({
    {{"default_num_steps", "1"},
     {"regrow_avb4", "true"}, {"particle_type", "0"}, {"site", "1"}, {"neighbor_index", "0"}, {"target_particle_type", "0"}, {"target_site", "1"}},
    {{"bond", "true"}, {"mobile_site", "0"}, {"anchor_site", "1"}},
    {{"angle", "true"}, {"mobile_site", "2"}, {"anchor_site", "0"}, {"anchor_site2", "1"}}}));
  mc.add(MakeTrialAdd({{"particle_type", "0"}}));
  mc.run(MakeRun({{"until_num_particles", "10"}}));
  mc.run(MakeRemove({{"name", "TrialAdd"}}));
  const std::string trials_per = feasst::str(1e5);
  mc.add(MakeEnergy());
//  mc.add(MakeLogAndMovie({{"output_file", "tmp/trimer2d"}, {"trials_per_write", trials_per}}));
  mc.add(MakeChirality2D());
  mc.add(MakeAnalyzeBonds());
  mc.add(MakeCheckEnergy({{"trials_per_update", trials_per}}));
  mc.add(MakeTune());
  auto mc2 = test_serialize_unique(mc);
  mc2->attempt(1e6);
  const Analyze& chiral = SeekAnalyze().reference("Chirality2D", *mc2);
  EXPECT_NEAR(chiral.accumulator().average(), 10., NEAR_ZERO);
  EXPECT_NEAR(chiral.accumulator().stdev(), 0., NEAR_ZERO);
  return mc2;
}

const double z_factor = 20.;

TEST(MonteCarlo, heterotrimer2d_VERY_LONG) {
  auto mc_no_avb = test_avb(false, false);
  Accumulator en_no_avb = SeekAnalyze().reference("Energy", *mc_no_avb).accumulator();
  INFO(en_no_avb.str());
  auto mc_avb2 = test_avb(true, false);
  Accumulator en_avb2 = SeekAnalyze().reference("Energy", *mc_avb2).accumulator();
  INFO(en_avb2.str());

  EXPECT_TRUE(en_no_avb.is_equivalent(en_avb2, z_factor, true));

  auto mc_avb4 = test_avb(true, false);
  Accumulator en_avb4 = SeekAnalyze().reference("Energy", *mc_avb4).accumulator();
  INFO(en_avb4.str());
  EXPECT_TRUE(en_no_avb.is_equivalent(en_avb4, z_factor, true));
}

TEST(MonteCarlo, multisite_neighbors) {
  std::vector<std::shared_ptr<EnergyMap> > neighs;
  neighs.push_back(MakeEnergyMapNeighbor());
  neighs.push_back(MakeEnergyMapNeighborCriteria());
  for (auto neigh : neighs) {
    INFO(neigh->class_name());
    MonteCarlo mc;
    //mc.set(MakeRandomMT19937({{"seed", "123"}}));
    //mc.set(MakeRandomMT19937({{"seed", "1757526510"}}));
    //mc.set(MakeRandomMT19937({{"seed", "1757524876"}}));
    mc.add(MakeConfiguration({{"cubic_side_length", "6"},
                              {"particle_type", "dimer:../particle/dimer.txt"}}));
    //mc.add(MakePotential(MakeLennardJones()));
    mc.add(MakeNeighborCriteria({{"maximum_distance", "1.5"}, {"minimum_distance", "1"}}));
    mc.add(MakePotential(MakeLennardJones(), MakeVisitModel(MakeVisitModelInner(neigh))));
    mc.add_to_reference(MakePotential(MakeLennardJones()));
    mc.set(MakeThermoParams({{"beta", "1"}, {"chemical_potential", "1"}}));
    mc.set(MakeMetropolis());
    //INFO("change trials per update?");
    mc.add(MakeCheckEnergy({{"trials_per_update", "1"}, {"decimal_places", "6"}}));
    mc.add(MakeTrialAdd({{"particle_type", "dimer"}}));
    //mc.add(MakeTrialTranslate({{"weight", "100"}})); // HWH remove?
    mc.run(MakeRun({{"until_num_particles", "5"}}));
    mc.run(MakeRemove({{"name", "TrialAdd"}}));
  //  auto neigh = MakeEnergyMapNeighbor();
  //  mc.set(0, MakePotential(MakeLennardJones(), MakeVisitModel(MakeVisitModelInner(neigh))));
  //  mc.add_to_reference(MakePotential(MakeLennardJones()));
    mc.set(MakeMetropolis());
  //  mc.add(MakeTrialTranslate());
  //  mc.add(MakeTrialRotate({{"tunable_param", "50"}}));
    mc.add(MakeTrialGrow({
      {{"default_num_steps", "4"}, {"default_reference_index", "0"},
       {"regrow", "true"}, {"particle_type", "dimer"}, {"site", "D1"}},
      {{"bond", "true"}, {"mobile_site", "D2"}, {"anchor_site", "D1"}}}));
    mc.add(MakeTrialGrow({
      {{"default_num_steps", "4"}, {"default_reference_index", "0"},
       {"regrow", "true"}, {"particle_type", "dimer"}, {"site", "D2"}},
      {{"bond", "true"}, {"mobile_site", "D1"}, {"anchor_site", "D2"}}}));
    EXPECT_EQ(4, mc.trial(0).stage(0).num_steps());
  //  mc.add(MakeLogAndMovie({{"trials_per_write", "100"}, {"output_file", "tmp/dimer"}}));
    for (int i = 0; i < 1e1; ++i) {
      mc.attempt(1);
      neigh->check(mc.configuration());
    }
    if (mc.system().potential(0).visit_model().inner().is_energy_map_queryable()) {
      INFO("Check:" << neigh->class_name());
      EXPECT_NEAR(mc.criteria().current_energy(), neigh->total_energy(), 1e-12);
    }
  }
}

void add_cg4_potential(MonteCarlo * mc, double eps_fc, double eps_fab) {
  Configuration * config = mc->get_system()->get_configuration();
  config->set_model_param("epsilon", 1, eps_fc);
  config->set_model_param("epsilon", 2, eps_fab);
  config->set_model_param("epsilon", 3, eps_fab);
  mc->add(MakePotential(MakeSquareWell()));
  { // intra is HardSphere with 90% reduced sigma
//    ModelParams params = mc->configuration().model_params().deep_copy();
    ModelParams params = mc->system().configuration().model_params().deep_copy();
    for (int type = 0; type < params.size(); ++type) {
      params.set("sigma", type, 0.9*params.select("sigma").value(type));
    }
    auto pot = MakePotential(MakeHardSphere(), MakeVisitModelIntra());
    pot->set(params);
    mc->add(pot);
  }
  // reference is HardSphere on all sites.
  mc->add_to_reference(MakePotential(MakeHardSphere()));
//  // reference is HardSphere on center with diameter 10
//  ModelParams params = mc->configuration().model_params();
//  for (const std::string parm : {"cutoff", "sigma"}) {
//    for (const int center : {0, 4}) params.set(parm, center, 10.);
//    for (const int branch : {1, 2, 3, 5, 6, 7}) params.set(parm, branch, 0.);
//  }
//  INFO(params.epsilon().str());
//  INFO(params.sigma().str());
//  INFO(params.cutoff().str());
//  auto ref = MakePotential(MakeHardSphere());
//  ref->set(params);
//  mc->add_to_reference(ref);
}

TEST(MonteCarlo, cg4_flexible_LONG) {
  MonteCarlo mc;
  mc.add(MakeConfiguration({
    {"cubic_side_length", "30"},
    {"particle_type0", "../plugin/chain/particle/cg4_mab_flex.txt"},
    {"add_particles_of_type0", "1"},
  }));
  EXPECT_EQ(1, mc.configuration().num_particles());
  EXPECT_EQ(1, mc.configuration().num_particles_of_type(0));
  add_cg4_potential(&mc, 1, 1);
  const double temperature = 0.7092;
  mc.set(MakeThermoParams({{"beta", str(1./temperature)}}));
  mc.set(MakeMetropolis());
//  mc.add(MakeTrialTranslate({{"reference_index", "0"}, {"tunable_param", "1"}}));
//  mc.add(MakeTrialRotate({{"reference_index", "0"}, {"tunable_param", "40"}}));
  //mc.add(MakeTrialGrow({{{"particle_type", "0"}, {"bond", "1"}, {"mobile_site", "1"}, {"anchor_site", "0"}}}));
  mc.add(MakeTrialGrow({{{"particle_type", "0"}, {"bond", "1"}, {"mobile_site", "1"}, {"anchor_site", "0"}, {"potential_acceptance", "1"}}}));
  mc.add(MakeTrialGrow({{{"particle_type", "0"}, {"bond", "1"}, {"mobile_site", "2"}, {"anchor_site", "0"}}}));
  mc.add(MakeTrialGrow({{{"particle_type", "0"}, {"bond", "1"}, {"mobile_site", "3"}, {"anchor_site", "0"}}}));
  std::string trials_per = "1e4";
  mc.add(MakeLog({{"trials_per_write", trials_per}, {"output_file", "tmp/cg4.txt"}}));
  mc.add(MakeMovie({{"trials_per_write", trials_per}, {"output_file", "tmp/cg4.xyz"}}));
  mc.add(MakeTune());
  auto bonds = MakeAnalyzeBonds({{"bond_bin_width", "0.05"}});
  mc.add(bonds);
  mc.attempt(1e6);

  // check that bonds extend to 7, and bonded particles do not overlap within 0.9sigma_ij
  //for (int type = 1; type < 2; ++type) {
  for (int type = 1; type < 4; ++type) {
    //INFO("type " << type);
    //INFO(bonds->bond(0).str());
    //INFO(bonds->bond_hist(type-1).str());
    EXPECT_NEAR(bonds->bond_hist(type-1).max(), 7.025, 1e-13);
    const double sigij = mc.system().potential(1).model_params().select("sigma").mixed_value(0, type);
    //INFO("sigij " << sigij);
    if (type == 1) {
      EXPECT_NEAR(0.9*3.91815, sigij, 0.00001);
    }
    int bin = bonds->bond_hist(type-1).bin(sigij);
    EXPECT_EQ(0, bonds->bond_hist(type-1).histogram()[bin - 1]);
    EXPECT_NE(0, bonds->bond_hist(type-1).histogram()[bin + 1]);
    bin = bonds->bond_hist(type-1).size() - 1;
    EXPECT_NEAR(0.5, bonds->bond_hist(type-1).histogram()[bin]/
                     bonds->bond_hist(type-1).histogram()[bin - 1], 0.1);
  }
}

void test_b2_cg4_flex(const bool mayer_intra) {
  MonteCarlo mc;
  //mc.set(MakeRandomMT19937({{"seed", "1663861442"}}));
  //mc.set(MakeRandomMT19937({{"seed", "1629905961"}}));
  //mc.set(MakeRandomMT19937({{"seed", "123"}}));
  auto config = MakeConfiguration({{"cubic_side_length", str(NEAR_INFINITY)},
    {"particle_type0", "../plugin/chain/particle/cg4_mab_flex.txt"}});
  config->add_particle_type("../plugin/chain/particle/cg4_mab_flex.txt");
  config->add_particle_of_type(0);
  config->add_particle_of_type(1);
  mc.add(config);
  EXPECT_EQ(2, mc.configuration().num_particles());
  EXPECT_EQ(1, mc.configuration().num_particles_of_type(0));
  add_cg4_potential(&mc, 1, 1);
  //const double temperature = 0.4688;
  const double temperature = 0.5309;
  //const double temperature = 0.7092;
  //add_cg4_potential(&mc, 1, 0);
  //const double temperature = 0.2097;
  //const double temperature = 0.2218;
  //const double temperature = 0.2438;
  //const double temperature = 0.2933;
  //add_cg4_potential(&mc, 1, 1.5);
  //const double temperature = 0.7090;
  //add_cg4_potential(&mc, 1, 2);
  //const double temperature = 0.8929;
  //add_cg4_potential(&mc, 0, 1);
  //const double temperature = 0.3769;
  mc.set(MakeThermoParams({{"beta", str(1./temperature)}}));
  auto mayer = MakeMayerSampling();
  int pot_acc = 1;
  if (mayer_intra) {
    mayer = MakeMayerSampling({{"intra_potential", "1"}});
    pot_acc = -1;
  }
  mc.set(mayer);
  for (const std::string ptype : {"0", "1"}) {
    std::string param = "0";
    if (ptype == "1") param = "1";
    mc.add(MakeTrialGrow({
      {{"default_reference_index", "0"}, {"default_new_only", "true"},
       {"particle_type", ptype}, {"translate", "true"}, {"site", "0"}, {"tunable_param", param}},
      {{"bond", "true"}, {"mobile_site", "1"}, {"anchor_site", "0"}, {"potential_acceptance", str(pot_acc)}},
      {{"bond", "true"}, {"mobile_site", "2"}, {"anchor_site", "0"}, {"potential_acceptance", str(pot_acc)}},
      {{"bond", "true"}, {"mobile_site", "3"}, {"anchor_site", "0"}, {"potential_acceptance", str(pot_acc)}}}));
  }
  std::string trials_per = "1e4";
  mc.add(MakeLog({{"trials_per_write", trials_per}, {"output_file", "tmp/cg4.txt"}}));
  mc.add(MakeMovie({{"trials_per_write", trials_per}, {"output_file", "tmp/cg4.xyz"}}));
  //mc.add(MakeCheckEnergy({{"trials_per", trials_per}}));
  //mc.add(MakeTune({{"trials_per", trials_per}}));
  auto bonds = MakeAnalyzeBonds({{"bond_bin_width", "0.05"}});
  mc.add(bonds);
  mc.attempt(1e6);
  INFO("mayer: " << mayer->mayer().str());
  INFO("mayer_ref: " << mayer->mayer_ref().str());
  INFO("b22 " << mayer->second_virial_ratio()
    << " +/- " << mayer->second_virial_ratio_block_stdev());
  EXPECT_NEAR(0, mayer->second_virial_ratio(), 8*mayer->second_virial_ratio_block_stdev());
//  INFO(bonds->bond_hist(0).str());
}

// HWH test that distributions are also unchanged whether or not potential_acceptance is used
TEST(MayerSampling, b2_cg4_flexible_LONG) {
  for (const bool mayer_intra : {true, false}) {
    INFO("mayer_intra " << mayer_intra);
    test_b2_cg4_flex(mayer_intra);
  }
}

// https://dx.doi.org/10.1063/1.4918557
TEST(MayerSampling, trimer_grow_LONG) {
  MonteCarlo mc;
  //mc.set(MakeRandomMT19937({{"seed", "123"}}));
  { auto config = MakeConfiguration({{"cubic_side_length", str(NEAR_INFINITY)}});
    config->add_particle_type("../particle/trimer_0.4L.txt");
//    config->add_particle_type("../particle/trimer_0.4L", "2.txt");
    config->add_particle_of_type(0);
    config->add_particle_of_type(0);
    //config->add_particle_of_type(1);
    const double rwca = std::pow(2, 1./6.);
    config->set_model_param("cutoff", 0, 1, rwca);
    config->set_model_param("cutoff", 1, 1, rwca);
//    config->set_model_param("cutoff", 0, 3, rwca);
//    config->set_model_param("cutoff", 1, 2, rwca);
//    config->set_model_param("cutoff", 2, 3, rwca);
    mc.add(config);
  }
  mc.add(MakePotential(MakeLennardJonesForceShift()));
  auto ref = MakePotential(MakeHardSphere());
  //ModelParams params = mc.system().configuration().model_params().deep_copy();
  auto params = ref->model_params(mc.system().configuration()).deep_copy();
  //params.set("sigma", 1, 0);
  params.set("sigma", 0, 1, 0);
  params.set("sigma", 1, 1, 0);
  ref->set(params);
  mc.add_to_reference(ref);
  mc.set(MakeThermoParams({{"beta", str(1./0.815)}}));
  auto mayer = MakeMayerSampling();
  mc.set(mayer);
//  mc.add(MakeTrialTranslate({{"new_only", "true"}, {"reference_index", "0"},
//    {"tunable_param", "1"}, {"particle_type", "0"}}));
//    //{"tunable_param", "1"}, {"particle_type", "1"}}));
//  mc.add(MakeTrialRotate({{"new_only", "true"}, {"reference_index", "0"},
//    {"tunable_param", "40"}}));
  for (const std::string ptype : {"0"}) {
  //for (const std::string ptype : {"0", "1"}) {
    mc.add(MakeTrialGrow({
      {{"default_reference_index", "0"}, {"default_new_only", "true"},
       {"particle_type", ptype}, {"translate", "true"}, {"site", "0"}},
      {{"bond", "true"}, {"mobile_site", "1"}, {"anchor_site", "0"}},
      {{"angle", "true"}, {"mobile_site", "2"}, {"anchor_site", "0"}, {"anchor_site2", "1"}}}));
//    mc.add(MakeTrialGrow({
//      {{"particle_type", ptype}, {"bond", "0"}, {"mobile_site", "0"}, {"anchor_site", "1"}},
//      {{"angle", "0"}, {"mobile_site", "2"}, {"anchor_site", "0"}, {"anchor_site2", "1"}},
//    }, {{"reference_index", "0"}, {"new_only", "true"}}));
  }
  const std::string trials_per = "1e4";
//  mc.add(MakeLogAndMovie({{"trials_per_write", trials_per}, {"output_file", "tmp/trib"}}));
  mc.add(MakeCheckEnergy({{"trials_per_update", trials_per}, {"tolerance", "1e-4"}}));
  mc.attempt(1e6);
  double b2hs = 2./3.*PI*std::pow(mc.configuration().model_params().select("sigma").value(0), 3); // A^3
  INFO(b2hs*mayer->second_virial_ratio());
  INFO("mayer: " << mayer->mayer().str());
  INFO("mayer_ref: " << mayer->mayer_ref().str());
  INFO("b22 " << mayer->second_virial_ratio()
    << " +/- " << mayer->second_virial_ratio_block_stdev());
  EXPECT_NEAR(0, mayer->mayer().average(), 4*mayer->mayer().block_stdev());
}

TEST(MonteCarlo, RigidBondAngleDihedral) {
  for (const std::string data : {
    "../particle/dimer.txt",
    "../particle/trimer_0.4L.txt",
    "../plugin/chain/test/data/tetramer_rigid.txt",
    }) {
    INFO(data);
    System system;
    system.add(MakeConfiguration({
      {"cubic_side_length", "10"},
      {"particle_type", "fluid:"+data},
      {"add_num_fluid_particles", "1"}}));
    system.set(MakeThermoParams({{"beta", "1"}}));
    system.precompute();
    auto random = MakeRandomMT19937();
    BondVisitor vis;

    std::shared_ptr<TrialSelect> select;
    std::shared_ptr<Perturb> perturb;
    if (data == "../particle/dimer.txt") {
      select = MakeTrialSelectBond({{"particle_type", "fluid"}, {"mobile_site", "D2"},
        {"anchor_site", "D1"}});
      perturb = MakePerturbDistance();
    } else if (data == "../particle/trimer_0.4L.txt") {
      select = MakeTrialSelectAngle({{"particle_type", "fluid"}, {"mobile_site", "2"},
        {"anchor_site", "0"}, {"anchor_site2", "1"}});
      perturb = MakePerturbDistanceAngle();
    } else if (data == "../plugin/chain/test/data/tetramer_rigid.txt") {
      select = MakeTrialSelectDihedral({{"particle_type", "fluid"}, {"mobile_site", "3"},
        {"anchor_site", "2"}, {"anchor_site2", "1"}, {"anchor_site3", "0"}});
      perturb = MakePerturbDihedral();
    } else {
      FATAL("unrecognized " << data);
    }
    select->precompute(&system);
    select->sel(&system, random.get());
    perturb->precompute(select.get(), &system);
    perturb->perturb(&system, select.get(), random.get());
    perturb->finalize(&system);
    FileXYZ().write_for_vmd("tmp/rigid.xyz", system.configuration());
    vis.compute_all(system.configuration());
    EXPECT_EQ(0, vis.energy());
  }
}

void test_equipartition(const std::string& data) {
  //for (const std::string num_steps : {"1"}) {
  for (const std::string num_steps : {"1", "4"}) {
    std::string ref = "-1";
    if (num_steps == "4") ref = "0";
    INFO("data " << data << " num_steps " << num_steps << " ref " << ref);
    MonteCarlo mc;
    //mc.set(MakeRandomMT19937({{"seed", "123"}}));
    mc.add(MakeConfiguration({
      {"particle_type0", "../plugin/chain/test/data/" + data + ".txt"},
      {"add_particles_of_type0", "1"},
      {"cubic_side_length", "10"}}));
    mc.add(MakePotential(MakeIdealGas()));
    mc.add_to_reference(MakePotential(MakeDontVisitModel()));
    mc.set(MakeThermoParams({{"beta", "1"}}));
    mc.set(MakeMetropolis());
    DEBUG("initial energy " << mc.criteria().current_energy());
    if (data == "dimer_harmonic") {
      mc.add(MakeTrialGrow({
        {{"default_num_steps", num_steps}, {"default_reference_index", ref},
         {"particle_type", "0"}, {"bond", "1"}, {"mobile_site", "1"}, {"anchor_site", "0"}}}));
    } else if (data == "trimer_rigid_angle" || data == "trimer_harmonic") {
      mc.add(MakeTrialGrow({
        {{"default_num_steps", num_steps}, {"default_reference_index", ref},
         {"particle_type", "0"}, {"bond", "1"}, {"mobile_site", "1"}, {"anchor_site", "0"}},
        {{"angle", "1"}, {"mobile_site", "2"}, {"anchor_site", "1"}, {"anchor_site2", "0"}}}));
    } else if (data == "tetramer_harmonic_no_dihedral") {
      mc.add(MakeTrialGrow({
        {{"default_num_steps", num_steps}, {"default_reference_index", ref},
         {"particle_type", "0"}, {"bond", "1"}, {"mobile_site", "1"}, {"anchor_site", "0"}},
        {{"angle", "1"}, {"mobile_site", "2"}, {"anchor_site", "1"}, {"anchor_site2", "0"}},
        {{"angle", "1"}, {"mobile_site", "3"}, {"anchor_site", "2"}, {"anchor_site2", "1"}}}));
    } else if (data == "tetramer_harmonic" ||
               data == "tetramer_harmonic_rigid_bond_angle" ||
               data == "tetramer_rigid") {
      mc.add(MakeTrialGrow({
        {{"default_num_steps", num_steps}, {"default_reference_index", ref},
         {"particle_type", "0"}, {"bond", "1"}, {"mobile_site", "1"}, {"anchor_site", "0"}},
        {{"angle", "1"}, {"mobile_site", "2"}, {"anchor_site", "1"}, {"anchor_site2", "0"}},
        {{"dihedral", "1"}, {"mobile_site", "3"}, {"anchor_site", "2"}, {"anchor_site2", "1"}, {"anchor_site3", "0"}}}));
    } else if (data == "tetramer_branched") {
      mc.add(MakeTrialGrow({
        {{"default_num_steps", num_steps}, {"default_reference_index", ref},
         {"particle_type", "0"}, {"bond", "1"}, {"mobile_site", "1"}, {"anchor_site", "0"}},
        {{"branch", "1"}, {"mobile_site", "2"}, {"mobile_site2", "3"}, {"anchor_site", "0"}, {"anchor_site2", "1"}}}));
    } else if (data == "pentamer_harmonic") {
      mc.add(MakeTrialGrow({
        {{"default_num_steps", num_steps}, {"default_reference_index", ref},
         {"particle_type", "0"}, {"bond", "1"}, {"mobile_site", "1"}, {"anchor_site", "0"}},
        {{"angle", "1"}, {"mobile_site", "2"}, {"anchor_site", "1"}, {"anchor_site2", "0"}},
        {{"dihedral", "1"}, {"mobile_site", "3"}, {"anchor_site", "2"}, {"anchor_site2", "1"}, {"anchor_site3", "0"}},
        {{"dihedral", "1"}, {"mobile_site", "4"}, {"anchor_site", "3"}, {"anchor_site2", "2"}, {"anchor_site3", "1"}}}));
    } else {
      FATAL("unrecognized: " << data);
    }
    //const std::string trials_per = "1";
    const std::string trials_per = "1e3";
//      mc.add(MakeLogAndMovie({{"trials_per_write", trials_per}, {"output_file", "tmp/harmonic"}}));
    mc.add(MakeCheckEnergy({{"trials_per_update", trials_per}}));
    auto en = MakeEnergy({{"trials_per_write", trials_per}});
    mc.add(en);
    auto bonds = MakeAnalyzeBonds({{"bond_bin_width", "0.05"}});
    // auto bonds = MakeAnalyzeBonds({{"bond_bin_width", "0.05"}, {"trials_per", "1e3"}});
    mc.add(bonds);
    mc.attempt(1e4);
    //INFO(bonds->bond_hist(0).str());
    //INFO(bonds->bond(0).average() << " +/- " << 3*bonds->bond(0).block_stdev());
    DEBUG(bonds->bond(0).str());
    double z_fac = 30;
    if (num_steps != "1") {
      z_fac *= 3;
      if (data == "tetramer_branched") { z_fac *= 2; }
    }
    double l_expect = 1. + 1/600.;  // 1 + 1/k
    if (data == "tetramer_harmonic_rigid_bond_angle" || data == "tetramer_rigid") {
      l_expect = 1.;
    }
    EXPECT_NEAR(l_expect, bonds->bond(0).average(), z_fac*bonds->bond(0).block_stdev());
    DEBUG(en->accumulator().str());
    double en_expect;
    if (data == "dimer_harmonic") {
      en_expect = 0.5;
    } else if (data == "trimer_rigid_angle") {
      en_expect = 1.;
    } else if (data == "trimer_harmonic") {
      en_expect = 1.5;
    } else if (data == "tetramer_harmonic_rigid_bond_angle") {
      en_expect = 0.5;
    } else if (data == "tetramer_harmonic_no_dihedral") {
      en_expect = 2.5;
    } else if (data == "tetramer_harmonic") {
      en_expect = 3.0;
    } else if (data == "tetramer_rigid") {
      en_expect = 0.0;
    } else if (data == "tetramer_branched") {
      en_expect = 3.0;
    } else if (data == "pentamer_harmonic") {
      en_expect = 4.5;
    } else {
      FATAL("unrecognized data: " << data);
    }
    EXPECT_NEAR(en_expect, en->accumulator().average(), z_fac*en->accumulator().block_stdev());
  }
}

// split up equipartition among different tests so it can run in parallel
TEST(MonteCarlo, equipartition_dimer_harmonic_LONG) { test_equipartition("dimer_harmonic"); }
TEST(MonteCarlo, equipartition_trimer_rigid_angle_LONG) { test_equipartition("trimer_rigid_angle"); }
TEST(MonteCarlo, equipartition_trimer_harmonic_LONG) { test_equipartition("trimer_harmonic"); }
TEST(MonteCarlo, equipartition_tetramer_harmonic_no_dihedral_LONG) { test_equipartition("tetramer_harmonic_no_dihedral"); }
TEST(MonteCarlo, equipartition_tetramer_harmonic_rigid_bond_angle_LONG) { test_equipartition("tetramer_harmonic_rigid_bond_angle"); }
TEST(MonteCarlo, equipartition_tetramer_rigid_LONG) { test_equipartition("tetramer_rigid"); }
TEST(MonteCarlo, equipartition_tetramer_harmonic_LONG) { test_equipartition("tetramer_harmonic"); }
TEST(MonteCarlo, equipartition_tetramer_branched_LONG) { test_equipartition("tetramer_branched"); }
TEST(MonteCarlo, equipartition_pentamer_harmonic_LONG) { test_equipartition("pentamer_harmonic"); }

TEST(MonteCarlo, single_butane) {
  MonteCarlo mc;
  //mc.set(MakeRandomMT19937({{"seed", "123"}}));
  mc.add(MakeConfiguration({
    {"particle_type0", "../particle/n-butane.txt"},
    {"add_particles_of_type0", "1"},
    {"cubic_side_length", "100"}}));
  mc.add(MakePotential(MakeLennardJones(),
                       MakeVisitModelIntra({{"intra_cut",  "3"}})));
  mc.set(MakeThermoParams({{"beta", "1"}}));
  mc.set(MakeMetropolis());
  DEBUG("initial energy " << mc.criteria().current_energy());
  mc.add(MakeTrialGrow({
    {{"particle_type", "0"}, {"bond", "1"}, {"mobile_site", "1"}, {"anchor_site", "0"}},
    {{"angle", "1"}, {"mobile_site", "2"}, {"anchor_site", "1"}, {"anchor_site2", "0"}},
    {{"dihedral", "1"}, {"mobile_site", "3"}, {"anchor_site", "2"}, {"anchor_site2", "1"}, {"anchor_site3", "0"}}}));
  //const std::string trials_per = "1";
  const std::string trials_per = "1e3";
//  mc.add(MakeLogAndMovie({{"trials_per_write", trials_per}, {"output_file", "tmp/butane"}}));
  mc.add(MakeCheckEnergy({{"trials_per_update", trials_per}}));
  auto en = MakeEnergy({{"trials_per_write", trials_per}, {"output_file", "tmp/but.csv"}});
  mc.add(en);
  const int bins = 20;
  auto bonds = MakeAnalyzeBonds({
      {"bond_bin_width", "0.05"}, {"angle_bin_width", "0.05"},
      {"dihedral_bin_width", str(PI/bins)}, {"dihedral_bin_center", str(PI/bins/2)}, {"output_file", "tmp/butbond.csv"}});
  // auto bonds = MakeAnalyzeBonds({{"bond_bin_width", "0.05"}, {"trials_per", "1e3"}});
  mc.add(bonds);
  mc.attempt(1e3);
  //INFO(bonds->bond_hist(0).str());
  //INFO(bonds->bond(0).average() << " +/- " << 3*bonds->bond(0).block_stdev());
//  INFO(bonds->bond(0).str());
//  INFO(bonds->bond_hist(0).str());
//  INFO(bonds->angle(0).str());
//  INFO(bonds->angle_hist(0).str());
//  INFO(bonds->dihedral(0).str());
//  INFO(bonds->dihedral_hist(0).str());
//  std::ofstream ss("tmp/butane_dihedral.txt");
//  ss << bonds->dihedral_hist(0).str();
  EXPECT_NEAR(bonds->dihedral_hist(0).histogram()[19], 700, 70);
}

TEST(MonteCarlo, ethane) {
  auto mc = MakeMonteCarlo();
  //mc->set(MakeRandomMT19937({{"seed", "123"}}));
  mc->add(MakeConfiguration({{"cubic_side_length", "30"},
                             {"particle_type0", "../particle/ethane.txt"},
                             {"cutoff", "15"}}));
  mc->add(MakePotential(MakeLennardJones()));
  mc->add(MakePotential(MakeLongRangeCorrections()));
  mc->set(MakeThermoParams({{"beta", str(1./325)},
                            {"chemical_potential", "-1"}}));
  mc->set(MakeMetropolis());
  mc->add(MakeTrialGrow({
    {{"bond", "true"},
     {"mobile_site", "0"},
     {"anchor_site", "1"},
     {"weight", "1"},
     {"enable_tunable", "true"},
     {"tunable_param", "1"},
     {"particle_type", "0"}}}));
  mc->add(MakeTrialGrow({
    {{"bond", "true"},
     {"mobile_site", "1"},
     {"anchor_site", "0"},
     {"weight", "1"},
     {"enable_tunable", "true"},
     {"tunable_param", "0.5"},
     {"particle_type", "0"}}}));
  mc->add(MakeTrialAdd({{"weight", "4"}, {"particle_type", "0"}}));
  mc->run(MakeRun({{"until_num_particles", "10"}}));
  mc->run(MakeRemove({{"name", "TrialAdd"}}));
  const std::string trials_per = "1e0";
  mc->add(MakeWrapParticles({{"trials_per_update", trials_per}}));
//  mc->add(MakeLogAndMovie({{"trials_per_write", trials_per}, {"output_file", "tmp/ethane.txt"}}));
  mc->add(MakeCheckEnergy({{"trials_per_update", trials_per}, {"tolerance", str(1e-2)}}));
  mc->add(MakeTune());
  mc->attempt(1e1);
}

TEST(MonteCarlo, water) {
  auto mc = MakeMonteCarlo();
  mc->set(MakeRandomMT19937({{"seed", "123"}}));
  mc->add(MakeConfiguration({{"cubic_side_length", "20"},
                             {"particle_type0", "../particle/spce.txt"},
                             {"cutoff", "10"}}));
  mc->add(MakePotential(MakeLennardJones(),
                        MakeVisitModelCutoffOuter()));
  mc->add(MakePotential(MakeLongRangeCorrections()));
  mc->set(MakeThermoParams({{"beta", str(1./325)},
                            {"chemical_potential", "-1"}}));
  mc->set(MakeMetropolis());
  mc->add(MakeTrialGrow({
    {{"angle", "true"},
     {"mobile_site", "1"},
     {"anchor_site", "0"},
     {"anchor_site2", "2"},
     {"weight", "1"},
     {"enable_tunable", "true"},
     {"tunable_param", "5"},
     {"particle_type", "0"}}}));
  mc->add(MakeTrialGrow({
    {{"angle", "true"},
     {"mobile_site", "2"},
     {"anchor_site", "0"},
     {"anchor_site2", "1"},
     {"weight", "1"},
     {"enable_tunable", "true"},
     {"tunable_param", "5"},
     {"particle_type", "0"}}}));
  mc->add(MakeTrialParticlePivot({{"particle_type", "0"}}));
  mc->add(MakeTrialCrankshaftSmall({
    {"particle_type", "0"},
    {"site", "0"},
    {"anchor_site0", "1"},
    {"anchor_site1", "2"},
    {"weight", "1."},
    {"tunable_param", "250."}}));
  mc->add(MakeTrialAdd({{"weight", "4"}, {"particle_type", "0"}}));
  mc->run(MakeRun({{"until_num_particles", "10"}}));
  mc->run(MakeRemove({{"name", "TrialAdd"}}));
  const std::string trials_per = "1e0";
  mc->add(MakeWrapParticles({{"trials_per_update", trials_per}}));
//  mc->add(MakeLogAndMovie({{"trials_per_write", trials_per}, {"output_file", "tmp/water"}}));
  mc->add(MakeCheckEnergy({{"trials_per_update", trials_per}, {"tolerance", str(1e-2)}}));
  mc->add(MakeTune());
  mc->attempt(1e1);
}

TEST(MonteCarlo, chainarglist) {
  auto mc = MakeMonteCarlo({{
    {"Configuration", {{"cubic_side_length", "20"},
                       {"particle_type", "../plugin/chain/particle/chain20.txt"}}},
    {"Potential", {{"Model", "IdealGas"}}},
    {"ThermoParams", {{"beta", "0.1"}, {"chemical_potential", "10"}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"tunable_param", "0.2"},
                        {"tunable_target_acceptance", "0.2"}}},
    {"TrialAdd", {{"particle_type", "0"}}},
    {"TrialGrowFile", {{"grow_file", "../plugin/chain/test/data/dimer_grow_file.txt"}}},
    {"Log", {{"trials_per_write", str(1e2)}, {"output_file", "tmp/lj.txt"}}},
    {"Movie", {{"trials_per_write", str(1e2)}, {"output_file", "tmp/lj.xyz"}}},
    {"CheckEnergy", {{"trials_per_update", str(1e2)}, {"tolerance", "1e-8"}}},
    {"Tune", {{}}},
//    {"Run", {{"until_num_particles", "50"}}},
//    {"ThermoParams", {{"beta", "1.2"}}},
//    {"Remove", {{"name", "TrialAdd"}}},
//    {"Run", {{"num_trials", str(1e3)}}},
//    {"Remove", {{"name", "Tune"}}},
//    {"Run", {{"num_trials", str(1e3)}}},
//    {"WriteCheckpoint", {{}}},
  }}, true);
}

TEST(MonteCarlo, angle_square_well) {
  auto mc = MakeMonteCarlo();
  mc->set(MakeRandomMT19937({{"seed", "123"}}));
  mc->add(MakeConfiguration({{"cubic_side_length", "20"},
                             {"particle_type0", "../plugin/chain/particle/chain3.txt"},
                             {"cutoff", "10"},
                             {"add_particles_of_type0", "1"},
                             {"group0", "end"},
                             {"end_site_type", "1"}}));
  mc->add(MakePotential(MakeLennardJones(),
                        MakeVisitModelCutoffOuter()));
  mc->set(MakeThermoParams({{"beta", "1"}, {"chemical_potential", "-1"}}));
  mc->set(MakeMetropolis());
  mc->add(MakeTrialGrow({
    {{"angle", "true"},
     {"mobile_site", "2"},
     {"anchor_site", "1"},
     {"anchor_site2", "0"},
     {"weight", "1"},
     {"particle_type", "0"}}}));
  const std::string trials_per = "1e0";
  mc->add(MakeMovie({{"trials_per_write", trials_per}, {"output_file", "tmp/chain3"}, {"group", "end"}}));
  mc->add(MakeCheckEnergy({{"trials_per_update", trials_per}, {"tolerance", str(1e-2)}}));
  for (int i = 0; i < 1e4; ++i) {
    mc->attempt(1);
    //std::cout << mc->configuration().particle(0).site(2).position().str() << std::endl;
  }
}

TEST(MonteCarlo, lj_position_swap) {
  auto mc = MakeMonteCarlo();
  mc->set(MakeRandomMT19937({{"seed", "123"}}));
  mc->add(MakeConfiguration({{"cubic_side_length", "8"},
    {"particle_type", "../particle/lj.txt,../particle/lj.txt"}}));
  mc->add(MakePotential(MakeLennardJones()));
  mc->add(MakePotential(MakeLongRangeCorrections()));
  mc->set(MakeThermoParams({{"beta", str(1./1.5)},
                            {"chemical_potential", "-1,-1"}}));
  mc->set(MakeMetropolis());
  mc->add(MakeTrialGrow({
    {{"position_swap", "true"},
     {"mobile_site", "0"},
     {"mobile_site2", "0"},
     {"particle_type", "0"},
     {"particle_type2", "1"},
     {"weight", "4"}}}));
  mc->add(MakeTrialAdd({{"weight", "1"}, {"particle_type", "0"}}));
  mc->run(MakeRun({{"until_num_particles", "2"}}));
  mc->run(MakeRemove({{"name", "TrialAdd"}}));
  mc->add(MakeTrialAdd({{"weight", "1"}, {"particle_type", "1"}}));
  mc->run(MakeRun({{"until_num_particles", "4"}}));
  mc->run(MakeRemove({{"name", "TrialAdd"}}));
  const std::string trials_per = "1e0";
  mc->add(MakeMovie({{"trials_per_write", trials_per}, {"output_file", "tmp/lj.xyz"}}));
  mc->add(MakeCheckEnergy({{"trials_per_update", trials_per}, {"tolerance", str(1e-2)}}));
  mc->add(MakeTune());
  mc->attempt(1e1);
}

//TEST(TrialGrow, reptate) {
//  MonteCarlo mc;
//  mc.add(MakeConfiguration({{"cubic_side_length", "20"},
//    {"particle_type0", "../plugin/chain/particle/chain5.txt"},
//    {"add_particles_of_type0", "1"}}));
//  mc.add(MakePotential(MakeLennardJones(), MakeVisitModelIntra({{"intra_cut", "1"}})));
//  mc.set(MakeThermoParams({{"beta", "1"}}));
//  mc.set(MakeMetropolis());
//  mc.add(MakeTrialGrow({
//    {{"reptate", "true"}, {"mobile_site", "0"}, {"anchor_site", "1"}, {"particle_type", "0"}},
//    {{"reptate", "true"}, {"mobile_site", "1"}, {"anchor_site", "2"}},
//    {{"reptate", "true"}, {"mobile_site", "2"}, {"anchor_site", "3"}},
//    {{"reptate", "true"}, {"mobile_site", "3"}, {"anchor_site", "4"}},
//    {{"bond", "true"}, {"mobile_site", "4"}, {"anchor_site", "3"}}}));
//  Particle chain = mc.configuration().particle(0);
//  mc.add(MakeMovie({{"output_file", "tmp/reptate.xyz"}}));
//  mc.add(MakeCheckEnergy());
//  while (mc.trial(0).acceptance() <= 0) {
//    mc.attempt(1);
//  }
//  for (int site = 0; site < 4; ++site) {
//    EXPECT_TRUE(chain.site(site+1).position().is_equal(mc.configuration().particle(0).site(site).position(), NEAR_ZERO));
//  }
//  mc.attempt(10);
//}

TEST(MonteCarlo, chain_npt_ideal_gas) {
  const std::string tpc = "1e1";
  auto mc = MakeMonteCarlo({{
    {"Configuration", {{"cubic_side_length", str(std::pow(5e4, 1./3))},
      {"particle_type", "../particle/n-butane.txt"},
      {"cutoff", "0"}}},
    {"Potential", {{"VisitModel", "DontVisitModel"}}},
    //{"Potential", {{"Model", "IdealGas"}}},
    //{"Potential", {{"Model", "LennardJones"}}},
    {"RefPotential", {{"VisitModel", "DontVisitModel"}}},
    {"ThermoParams", {{"beta", "1"}, {"chemical_potential", "1."}, {"pressure", "1e-3"}}},
    {"Metropolis", {{}}},
    {"TrialGrowFile", {{"grow_file", "../plugin/chain/test/data/c4_grow_ce.txt"}}},
    {"TrialGrowFile", {{"grow_file", "../plugin/chain/test/data/c4_grow_gce.txt"}}},
    {"Remove", {{"name_contains", "remove"}}},
    {"Run", {{"until_num_particles", "50"}}},
    {"Remove", {{"name_contains", "add"}}},
    {"TrialVolume", {{"reference_index", "0"}, {"weight", "1e1"}, {"tunable_param", "0.8"}}},
    {"Metropolis", {{"trials_per_cycle", tpc}, {"cycles_to_complete", "1e2"}}},
    {"Log", {{"trials_per_write", tpc}, {"output_file", "tmp/lj.csv"}}},
    {"Movie", {{"trials_per_write", tpc}, {"output_file", "tmp/lj.xyz"}}},
    {"Volume", {{"trials_per_write", tpc}, {"output_file", "tmp/lj_vol.csv"}}},
    {"CheckEnergy", {{"trials_per_update", tpc}, {"decimal_places", "7"}}},
    {"Tune", {{"trials_per_tune", "20"}}},
    {"Run", {{"until", "complete"}}},
  }}, true);
  const Accumulator& vol = mc->analyzers().back()->accumulator();
  DEBUG("vol:" << vol.str());
  EXPECT_NEAR(5e4, vol.average(), 5*vol.block_stdev());
}

TEST(MonteCarlo, gibbs_trappe_ideal_gas_LONG) {
  const std::string tpc = "1e3";
  auto mc = MakeMonteCarlo({{
    {"Configuration", {{"cubic_side_length", "40"}, {"particle_type", "../particle/n-butane.txt"}, {"cutoff", "0"}}},
    {"Configuration", {{"cubic_side_length", "40"}, {"particle_type", "../particle/n-butane.txt"}, {"cutoff", "0"}}},
    {"CopyFollowingLines", {{"for_num_configurations", "2"}}},
    {"Potential", {{"Model", "IdealGas"}}},
    {"RefPotential", {{"VisitModel", "DontVisitModel"}}},
    {"EndCopy", {{}}},
    {"ThermoParams", {{"beta", "1"}, {"chemical_potential", "10"}}},
    {"Metropolis", {{}}},

    // fill box
    {"TrialGrowFile", {{"grow_file", "../plugin/chain/test/data/c4_grow_ce.txt"}}},
    {"TrialGrowFile", {{"grow_file", "../plugin/chain/test/data/c4_grow_gce.txt"}}},
    {"Remove", {{"name_contains", "remove"}}},
    {"Run", {{"until_num_particles", "50"}, {"configuration_index", "0"}}},
    {"Remove", {{"name_contains", "add"}}},
    {"TrialGrowFile", {{"grow_file", "../plugin/chain/test/data/c4_grow_gce_conf1.txt"}}},
    {"Remove", {{"name_contains", "remove"}}},
    {"Run", {{"until_num_particles", "50"}, {"configuration_index", "1"}}},
    {"Remove", {{"name_contains", "add"}}},

    {"Metropolis", {{"trials_per_cycle", tpc}, {"cycles_to_complete", "1e2"}}},
    {"TrialGrowFile", {{"grow_file", "../plugin/chain/test/data/c4_grow_gibbs.txt"}}},
    {"TrialGibbsVolumeTransfer", {{"tunable_param", "10"}, {"reference_index", "0"}}},
    {"Tune", {{"trials_per_tune", "4"}}},
    {"Log", {{"trials_per_write", tpc}, {"output_file", "tmp/lj.csv"}}},
    {"CopyFollowingLines", {{"for_num_configurations", "2"}, {"replace_with_index", "[config]"}}},
    {"Movie", {{"trials_per_write", tpc}, {"output_file", "tmp/ljc[config].xyz"}}},
    {"Density", {{"trials_per_write", tpc}, {"output_file", "tmp/ljc[config]dens.csv"}}},
    {"NumParticles", {{"trials_per_write", tpc}, {"output_file", "tmp/ljc[config]num.csv"}}},
    {"Volume", {{"trials_per_write", tpc}, {"output_file", "tmp/ljc[config]vol.csv"}}},
    {"EndCopy", {{}}},
    {"CheckEnergy", {{"trials_per_update", tpc}, {"tolerance", str(1e-9)}}},
    //{"GibbsInitialize", {{"trials_per_update", tpc}}},
    {"Run", {{"until", "complete"}}},
  //}});
  }}, true);
  auto mc2 = test_serialize_unique(*mc);
  EXPECT_EQ(2, mc2->system().num_configurations());
  const Accumulator& num0 = mc2->analyze(5).accumulator();
  const Accumulator& num1 = mc2->analyze(6).accumulator();
  DEBUG(num0.str() << std::endl << num1.str());
  EXPECT_NEAR(50, num0.average(), 5*num0.block_stdev());
  EXPECT_NEAR(50, num1.average(), 5*num1.block_stdev());
  const Accumulator& vol0 = mc2->analyze(7).accumulator();
  const Accumulator& vol1 = mc2->analyze(8).accumulator();
  DEBUG(vol0.str() << std::endl << vol1.str());
  EXPECT_NEAR(std::pow(40, 3), vol0.average(), 5*vol0.block_stdev());
  EXPECT_NEAR(std::pow(40, 3), vol1.average(), 5*vol1.block_stdev());
}

}  // namespace feasst
