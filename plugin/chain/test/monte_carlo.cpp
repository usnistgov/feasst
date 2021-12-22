#include <memory>
#include "utils/test/utils.h"
#include "utils/include/io.h"
#include "math/include/accumulator.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/domain.h"
#include "system/include/ideal_gas.h"
#include "system/include/hard_sphere.h"
#include "system/include/lennard_jones.h"
#include "system/include/long_range_corrections.h"
#include "system/include/visit_model_intra_map.h"
#include "system/include/dont_visit_model.h"
#include "models/include/square_well.h"
#include "models/include/lennard_jones_force_shift.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trials.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/run.h"
#include "monte_carlo/include/seek_num_particles.h"
#include "monte_carlo/include/trial_compute_move.h"
#include "monte_carlo/include/trial_select_dihedral.h"
#include "monte_carlo/include/perturb_dihedral.h"
#include "steppers/include/log.h"
#include "steppers/include/check_energy_and_tune.h"
#include "steppers/include/log_and_movie.h"
#include "steppers/include/movie.h"
#include "steppers/include/tune.h"
#include "steppers/include/check.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/check_physicality.h"
#include "steppers/include/check_properties.h"
#include "steppers/include/chirality_2d.h"
#include "steppers/include/energy.h"
#include "steppers/include/seek_analyze.h"
#include "cluster/include/energy_map_all.h"
#include "cluster/include/trial_avb2.h"
#include "cluster/include/energy_map_neighbor.h"
#include "chain/include/trials.h"
#include "chain/include/trial_grow.h"
#include "chain/include/trial_grow_linear.h"
#include "chain/include/analyze_bonds.h"
#include "chain/test/system_chain.h"
#include "mayer/include/mayer_sampling.h"

namespace feasst {

TEST(MonteCarlo, chain) {
  MonteCarlo mc;
  mc.set(MakeRandomMT19937({{"seed", "1610724059"}}));
  mc.set(chain_system());
  mc.set(MakeThermoParams({{"beta", "1"}, {"chemical_potential", "1."}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialAdd({{"particle_type", "0"}}));
  mc.run(MakeRun({{"until_num_particles", "1"}}));
  mc.run(MakeRemoveTrial({{"name", "TrialAdd"}}));
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
  const int steps_per = 1e0;
  mc.add(MakeLog({
    {"steps_per", str(steps_per)},
    {"file_name", "tmp/chainlog.txt"}}));
  mc.add(MakeMovie({
    {"steps_per", str(steps_per)},
    {"file_name", "tmp/chain10movie.xyz"}}));
  mc.add(MakeCheckEnergy({
    {"steps_per", str(steps_per)},
    {"tolerance", "1e-10"}}));
  mc.add(MakeTune({{"steps_per", str(steps_per)}}));
  mc.attempt(3e2);

  MonteCarlo mc2 = test_serialize(mc);
  EXPECT_EQ(mc2.analyzers().size(), 2);
}

// HWH this test is known to fail infrequently
// HWH sometimes it happens with valgrind and serialization, which is why the test is labelled long
TEST(MonteCarlo, TrialGrow_LONG) {
  for (const std::string particle : {"lj", "spce"}) {
    INFO(particle);
    double box_length = 8.;
    std::string data = "../forcefield/dimer.fstprt";
    if (particle == "spce") {
      box_length=20;
      data = "../forcefield/spce.fstprt";
    }
    MonteCarlo mc;
    //mc.set(MakeRandomMT19937({{"seed", "1635356012"}}));
    mc.add(MakeConfiguration({{"cubic_box_length", str(box_length)},
                              {"particle_type0", data}}));
    mc.add(MakePotential(MakeLennardJones()));
    mc.add(MakePotential(MakeLongRangeCorrections()));
    mc.add_to_reference(MakePotential(MakeLennardJones()));
    mc.set(MakeThermoParams({{"beta", "1.2"}, {"chemical_potential", "1."}}));
    mc.set(MakeMetropolis());
    mc.add(MakeTrialAdd({{"particle_type", "0"}}));
    mc.run(MakeRun({{"until_num_particles", "3"}}));
    mc.run(MakeRemoveTrial({{"name", "TrialAdd"}}));
    mc.set(MakeThermoParams({{"beta", "1.2"}, {"chemical_potential", "-700"}}));
    std::vector<argtype> grow_args = {
      {{"transfer", "true"},
       {"regrow", "true"},
       {"particle_type", "0"},
       {"site", "0"},
       {"weight", "100"}},
      {{"bond", "true"},
       {"mobile_site", "1"},
       {"anchor_site", "0"},
       {"num_steps", "5"}}};
    if (particle == "spce") grow_args.push_back(
      {{"angle", "true"},
       {"mobile_site", "2"},
       {"anchor_site", "0"},
       {"anchor_site2", "1"}});
    mc.add(MakeTrialGrow(grow_args, {{"num_steps", "4"}, {"reference_index", "0"}}));
    EXPECT_EQ(4, mc.trial(0).stage(0).rosenbluth().num());
    EXPECT_EQ(5, mc.trial(0).stage(1).rosenbluth().num());
    mc.add(MakeLogAndMovie({{"steps_per", str(1e0)}, {"file_name", "tmp/lj"}}));
    mc.add(MakeCheckEnergyAndTune({{"steps_per", str(1e0)}, {"tolerance", str(1e-9)}}));
    EXPECT_EQ(3, mc.trials().num());
    EXPECT_TRUE(mc.trial(0).stage(0).trial_select().is_ghost());   // add
    EXPECT_FALSE(mc.trial(1).stage(0).trial_select().is_ghost());  // remove
    EXPECT_FALSE(mc.trial(2).stage(0).trial_select().is_ghost());  // regrow
    mc.add(MakeMovie(
     {{"steps_per", "1"},
      {"file_name", "tmp/grow.xyz"}}));
    for (int trial = 0; trial < 2e1; ++trial) {
      mc.attempt(1);
      //mc.configuration().check();
    }
    EXPECT_LT(mc.configuration().num_particles(), 3);
    mc.set(MakeThermoParams({{"beta", "1.2"}, {"chemical_potential", "100"}}));
    MonteCarlo mc2 = test_serialize(mc);
    mc2.attempt(2e1);
    EXPECT_GE(mc2.configuration().num_particles(), 1);
    mc2.configuration().check();
    // INFO(mc.trial(1)->accept().perturbed().str());
  }
}

MonteCarlo cg7mab2(const std::string& data, const int num, const int steps_per = 1) {
  INFO("data " << data);
  MonteCarlo mc;
  //mc.set(MakeRandomMT19937({{"seed", "1633624429"}}));
  mc.add(MakeConfiguration({{"cubic_box_length", "30"},
    {"particle_type", "../plugin/chain/forcefield/" + data},
    {"set_cutoff_min_to_sigma", "true"}}));
  mc.add(MakePotential(MakeHardSphere()));
  if (is_found_in(data, "fullangflex")) {
    mc.add(MakePotential(MakeHardSphere(),
                     MakeVisitModelIntraMap({{"exclude_bonds", "true"}})));
  }
  mc.set(MakeThermoParams({{"beta", "1."}, {"chemical_potential", "1"}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialAdd({{"particle_type", "0"}}));
  mc.run(MakeRun({{"until_num_particles", str(num)}}));
  mc.run(MakeRemoveTrial({{"name", "TrialAdd"}}));
  if (is_found_in(data, "fullangflex")) {
    mc.add(MakeTrialGrow({
      //{{"particle_type", "0"}, {"site", "0"}, {"weight", "4"}, {"regrow", "1"}},
      //{{"bond", "1"}, {"mobile_site", "2"}, {"anchor_site", "0"}},
      {{"particle_type", "0"}, {"weight", "4"}, {"bond", "1"}, {"mobile_site", "1"}, {"anchor_site", "0"}},
      {{"angle", "1"}, {"mobile_site", "2"}, {"anchor_site", "1"}, {"anchor_site2", "0"}},
      {{"branch", "1"}, {"mobile_site", "3"}, {"mobile_site2", "5"}, {"anchor_site", "0"}, {"anchor_site2", "1"}},
      {{"angle", "1"}, {"mobile_site", "4"}, {"anchor_site", "3"}, {"anchor_site2", "0"}},
      {{"angle", "1"}, {"mobile_site", "6"}, {"anchor_site", "5"}, {"anchor_site2", "0"}}}));
  } else {
    mc.add(MakeTrialGrow({
      //{{"particle_type", "0"}, {"site", "0"}, {"weight", "4"}, {"regrow", "1"}},
      //{{"bond", "1"}, {"mobile_site", "2"}, {"anchor_site", "0"}},
      {{"particle_type", "0"}, {"weight", "4"}, {"bond", "1"}, {"mobile_site", "2"}, {"anchor_site", "0"}},
      {{"angle", "1"}, {"mobile_site", "1"}, {"anchor_site", "2"}, {"anchor_site2", "0"}},
      {{"branch", "1"}, {"mobile_site", "4"}, {"mobile_site2", "6"}, {"anchor_site", "0"}, {"anchor_site2", "2"}},
      {{"angle", "1"}, {"mobile_site", "3"}, {"anchor_site", "4"}, {"anchor_site2", "0"}},
      {{"angle", "1"}, {"mobile_site", "5"}, {"anchor_site", "6"}, {"anchor_site2", "0"}}}));
  }
  mc.add(MakeLogAndMovie({{"steps_per", str(steps_per)}, {"file_name", "tmp/" + data}}));
  mc.add(MakeCheckEnergyAndTune({{"steps_per", str(steps_per)}, {"tolerance", str(1e-9)}}));
  return mc;
}

TEST(MonteCarlo, cg7mab2) {
  for (const std::string data : {
      "cg7mab2.fstprt",
      "cg7mab2closed.fstprt",
      "cg7mab2extended.fstprt",
      "cg7mab2flex.fstprt",
      "cg7mab2fullangflex.fstprt",
    }) {
    cg7mab2(data, 1).attempt(1e2);
  }
}

TEST(MonteCarlo, cg7mab2_LONG) {
  for (const std::string data : {
      "cg7mab2.fstprt",
      "cg7mab2closed.fstprt",
      "cg7mab2extended.fstprt",
      "cg7mab2flex.fstprt",
      "cg7mab2fullangflex.fstprt",
    }) {
    cg7mab2(data, 10, 1e4).attempt(1e6);
  }
}

TEST(System, Angles2D) {
  MonteCarlo mc;
  mc.add(MakeConfiguration({{"side_length0", "6"}, {"side_length1", "6"},
    {"particle_type", "../plugin/chain/forcefield/heterotrimer2d.fstprt"},
    {"add_particles_of_type0", "1"}}));
  mc.add(MakePotential(MakeLennardJones()));
  mc.set(MakeThermoParams({{"beta", "1"}}));
  mc.set(MakeMetropolis());
  INFO(mc.criteria().current_energy());
}

MonteCarlo test_avb(const bool avb2, const bool avb4 = true) {
  MonteCarlo mc;
  mc.set(MakeRandomMT19937({{"seed", "time"}}));
  //mc.set(MakeRandomMT19937({{"seed", "123"}}));
  mc.add(MakeConfiguration({{"side_length0", "6"}, {"side_length1", "6"},
    {"particle_type", "../plugin/chain/forcefield/heterotrimer2d.fstprt"}}));
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
    {{"regrow_avb2", "true"}, {"particle_type", "0"}, {"site", "1"}, {"neighbor_index", "0"}, {"target_particle_type", "0"}, {"target_site", "1"}},
    {{"bond", "true"}, {"mobile_site", "0"}, {"anchor_site", "1"}},
    {{"angle", "true"}, {"mobile_site", "2"}, {"anchor_site", "0"}, {"anchor_site2", "1"}}},
    {{"num_steps", "1"}}));
  if (avb4) mc.add(MakeTrialGrow({
    {{"regrow_avb4", "true"}, {"particle_type", "0"}, {"site", "1"}, {"neighbor_index", "0"}, {"target_particle_type", "0"}, {"target_site", "1"}},
    {{"bond", "true"}, {"mobile_site", "0"}, {"anchor_site", "1"}},
    {{"angle", "true"}, {"mobile_site", "2"}, {"anchor_site", "0"}, {"anchor_site2", "1"}}},
    {{"num_steps", "1"}}));
  mc.add(MakeTrialAdd({{"particle_type", "0"}}));
  mc.run(MakeRun({{"until_num_particles", "10"}}));
  mc.run(MakeRemoveTrial({{"name", "TrialAdd"}}));
  const std::string steps_per = feasst::str(1e5);
  mc.add(MakeEnergy());
  mc.add(MakeLogAndMovie({{"file_name", "tmp/trimer2d"}, {"steps_per", steps_per}}));
  mc.add(MakeChirality2D());
  mc.add(MakeAnalyzeBonds());
  mc.add(MakeCheckEnergyAndTune({{"steps_per", steps_per}}));
  MonteCarlo mc2 = test_serialize(mc);
  mc2.attempt(1e6);
  const Analyze& chiral = SeekAnalyze().reference("Chirality2D", mc2);
  EXPECT_NEAR(chiral.accumulator().average(), 10., NEAR_ZERO);
  EXPECT_NEAR(chiral.accumulator().stdev(), 0., NEAR_ZERO);
  return mc2;
}

const double z_factor = 20.;

TEST(MonteCarlo, heterotrimer2d_VERY_LONG) {
  MonteCarlo mc_no_avb = test_avb(false, false);
  Accumulator en_no_avb = SeekAnalyze().reference("Energy", mc_no_avb).accumulator();
  INFO(en_no_avb.str());
  MonteCarlo mc_avb2 = test_avb(true, false);
  Accumulator en_avb2 = SeekAnalyze().reference("Energy", mc_avb2).accumulator();
  INFO(en_avb2.str());

  EXPECT_TRUE(en_no_avb.is_equivalent(en_avb2, z_factor, true));

  MonteCarlo mc_avb4 = test_avb(true, false);
  Accumulator en_avb4 = SeekAnalyze().reference("Energy", mc_avb4).accumulator();
  INFO(en_avb4.str());
  EXPECT_TRUE(en_no_avb.is_equivalent(en_avb4, z_factor, true));
}

TEST(MonteCarlo, multisite_neighbors) {
  MonteCarlo mc;
  //mc.set(MakeRandomMT19937({{"seed", "123"}}));
  //mc.set(MakeRandomMT19937({{"seed", "1610132694"}}));
  mc.add(MakeConfiguration({{"cubic_box_length", "6"},
                            {"particle_type0", "../forcefield/dimer.fstprt"}}));
  mc.add(MakePotential(MakeLennardJones()));
  mc.set(MakeThermoParams({{"beta", "1"}, {"chemical_potential", "1"}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialAdd({{"particle_type", "0"}}));
  mc.run(MakeRun({{"until_num_particles", "5"}}));
  mc.run(MakeRemoveTrial({{"name", "TrialAdd"}}));
  auto neigh = MakeEnergyMapNeighbor();
  mc.set(0, MakePotential(MakeLennardJones(), MakeVisitModel(MakeVisitModelInner(neigh))));
  mc.add_to_reference(MakePotential(MakeLennardJones()));
  mc.set(MakeMetropolis());
//  mc.add(MakeTrialTranslate());
//  mc.add(MakeTrialRotate({{"tunable_param", "50"}}));
  mc.add(MakeTrialGrow({
    {{"regrow", "true"}, {"particle_type", "0"}, {"site", "0"}},
    {{"bond", "true"}, {"mobile_site", "1"}, {"anchor_site", "0"}}},
    {{"num_steps", "4"}, {"reference_index", "0"}}));
  mc.add(MakeTrialGrow({
    {{"regrow", "true"}, {"particle_type", "0"}, {"site", "1"}},
    {{"bond", "true"}, {"mobile_site", "0"}, {"anchor_site", "1"}}},
    {{"num_steps", "4"}, {"reference_index", "0"}}));
  mc.add(MakeLogAndMovie({{"steps_per", "100"}, {"file_name", "tmp/dimer"}}));
  for (int i = 0; i < 1e1; ++i) {
    mc.attempt(1);
    neigh->check(mc.configuration());
  }
  EXPECT_NEAR(mc.criteria().current_energy(), neigh->total_energy(), 1e-12);
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
    {"cubic_box_length", "30"},
    {"particle_type0", install_dir() + "/plugin/chain/forcefield/cg4_mab_flex.fstprt"},
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
  std::string steps_per = "1e4";
  mc.add(MakeLog({{"steps_per", steps_per}, {"file_name", "tmp/cg4.txt"}}));
  mc.add(MakeMovie({{"steps_per", steps_per}, {"file_name", "tmp/cg4.xyz"}}));
  mc.add(MakeTune({{"steps_per", steps_per}}));
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

// HWH test that distributions are also unchanged whether or not potential_acceptance is used
TEST(MayerSampling, b2_cg4_flexible_LONG) {
  MonteCarlo mc;
  //mc.set(MakeRandomMT19937({{"seed", "1629905961"}}));
  //mc.set(MakeRandomMT19937({{"seed", "123"}}));
  auto config = MakeConfiguration({{"cubic_box_length", str(NEAR_INFINITY)},
    {"particle_type0", install_dir() + "/plugin/chain/forcefield/cg4_mab_flex.fstprt"}});
  config->add_particle_type(install_dir() + "/plugin/chain/forcefield/cg4_mab_flex.fstprt", "2");
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
  mc.set(mayer);
  for (const std::string ptype : {"0", "1"}) {
    std::string param = "0";
    if (ptype == "1") param = "1";
    mc.add(MakeTrialGrow({
      {{"particle_type", ptype}, {"translate", "true"}, {"site", "0"}, {"tunable_param", param}},
      {{"bond", "true"}, {"mobile_site", "1"}, {"anchor_site", "0"}, {"potential_acceptance", "1"}},
      {{"bond", "true"}, {"mobile_site", "2"}, {"anchor_site", "0"}, {"potential_acceptance", "1"}},
      {{"bond", "true"}, {"mobile_site", "3"}, {"anchor_site", "0"}, {"potential_acceptance", "1"}}
    }, {{"reference_index", "0"}, {"new_only", "true"}}));
  }
  std::string steps_per = "1e4";
  mc.add(MakeLog({{"steps_per", steps_per}, {"file_name", "tmp/cg4.txt"}}));
  mc.add(MakeMovie({{"steps_per", steps_per}, {"file_name", "tmp/cg4.xyz"}}));
  mc.add(MakeCheckEnergy({{"steps_per", steps_per}}));
  //mc.add(MakeTune({{"steps_per", steps_per}}));
  auto bonds = MakeAnalyzeBonds({{"bond_bin_width", "0.05"}});
  mc.add(bonds);
  mc.attempt(1e7);
  INFO("mayer: " << mayer->mayer().str());
  INFO("mayer_ref: " << mayer->mayer_ref().str());
  INFO("b22 " << mayer->second_virial_ratio()
    << " +/- " << mayer->second_virial_ratio_block_stdev());
  EXPECT_NEAR(0, mayer->second_virial_ratio(), 8*mayer->second_virial_ratio_block_stdev());
//  INFO(bonds->bond_hist(0).str());
}

// https://dx.doi.org/10.1063/1.4918557
TEST(MayerSampling, trimer_grow_LONG) {
  MonteCarlo mc;
  //mc.set(MakeRandomMT19937({{"seed", "123"}}));
  { auto config = MakeConfiguration({{"cubic_box_length", str(NEAR_INFINITY)}});
    config->add_particle_type(install_dir() + "/forcefield/trimer_0.4L.fstprt");
//    config->add_particle_type(install_dir() + "/forcefield/trimer_0.4L", "2.fstprt");
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
      {{"particle_type", ptype}, {"translate", "true"}, {"site", "0"}},
      {{"bond", "true"}, {"mobile_site", "1"}, {"anchor_site", "0"}},
      {{"angle", "true"}, {"mobile_site", "2"}, {"anchor_site", "0"}, {"anchor_site2", "1"}},
    }, {{"reference_index", "0"}, {"new_only", "true"}}));
//    mc.add(MakeTrialGrow({
//      {{"particle_type", ptype}, {"bond", "0"}, {"mobile_site", "0"}, {"anchor_site", "1"}},
//      {{"angle", "0"}, {"mobile_site", "2"}, {"anchor_site", "0"}, {"anchor_site2", "1"}},
//    }, {{"reference_index", "0"}, {"new_only", "true"}}));
  }
  const std::string steps_per = "1e4";
  mc.add(MakeLogAndMovie({{"steps_per", steps_per}, {"file_name", "tmp/trib"}}));
  mc.add(MakeCheckEnergy({{"steps_per", steps_per}, {"tolerance", "1e-4"}}));
  mc.attempt(1e6);
  double b2hs = 2./3.*PI*std::pow(mc.configuration().model_params().select("sigma").value(0), 3); // A^3
  INFO(b2hs*mayer->second_virial_ratio());
  INFO("mayer: " << mayer->mayer().str());
  INFO("mayer_ref: " << mayer->mayer_ref().str());
  INFO("b22 " << mayer->second_virial_ratio()
    << " +/- " << mayer->second_virial_ratio_block_stdev());
  EXPECT_NEAR(0, mayer->mayer().average(), 4*mayer->mayer().block_stdev());
}

TEST(TrialGrow, reptate) {
  MonteCarlo mc;
  mc.add(MakeConfiguration({{"cubic_box_length", "20"},
    {"particle_type0", "../plugin/chain/forcefield/chain5.fstprt"},
    {"add_particles_of_type0", "1"}}));
  mc.add(MakePotential(MakeLennardJones(), MakeVisitModelIntra({{"cutoff", "1"}})));
  mc.set(MakeThermoParams({{"beta", "1"}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialGrow({
    {{"reptate", "true"}, {"mobile_site", "0"}, {"anchor_site", "1"}, {"particle_type", "0"}},
    {{"reptate", "true"}, {"mobile_site", "1"}, {"anchor_site", "2"}},
    {{"reptate", "true"}, {"mobile_site", "2"}, {"anchor_site", "3"}},
    {{"reptate", "true"}, {"mobile_site", "3"}, {"anchor_site", "4"}},
    {{"bond", "true"}, {"mobile_site", "4"}, {"anchor_site", "3"}}}));
  Particle chain = mc.configuration().particle(0);
  mc.add(MakeMovie({{"file_name", "tmp/reptate.xyz"}}));
  mc.add(MakeCheckEnergy());
  while (mc.trial(0).acceptance() <= 0) {
    mc.attempt(1);
  }
  for (int site = 0; site < 4; ++site) {
    EXPECT_TRUE(chain.site(site+1).position().is_equal(mc.configuration().particle(0).site(site).position(), NEAR_ZERO));
  }
  mc.attempt(10);
}

TEST(MonteCarlo, RigidBondAngleDihedral) {
  for (const std::string data : {
    "../forcefield/dimer.fstprt",
    "../forcefield/trimer_0.4L.fstprt",
    "../plugin/chain/test/data/tetramer_rigid.fstprt",
    }) {
    INFO(data);
    System system;
    system.add(*MakeConfiguration({
      {"cubic_box_length", "10"},
      {"particle_type", data},
      {"add_particles_of_type0", "1"}}));
    system.set(MakeThermoParams({{"beta", "1"}}));
    system.precompute();
    auto random = MakeRandomMT19937();
    BondVisitor vis;

    std::shared_ptr<TrialSelect> select;
    std::shared_ptr<Perturb> perturb;
    if (data == "../forcefield/dimer.fstprt") {
      select = MakeTrialSelectBond({{"particle_type", "0"}, {"mobile_site", "1"},
        {"anchor_site", "0"}});
      perturb = MakePerturbDistance();
    } else if (data == "../forcefield/trimer_0.4L.fstprt") {
      select = MakeTrialSelectAngle({{"particle_type", "0"}, {"mobile_site", "2"},
        {"anchor_site", "0"}, {"anchor_site2", "1"}});
      perturb = MakePerturbDistanceAngle();
    } else if (data == "../plugin/chain/test/data/tetramer_rigid.fstprt") {
      select = MakeTrialSelectDihedral({{"particle_type", "0"}, {"mobile_site", "3"},
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

TEST(MonteCarlo, equipartition_LONG) {
  for (const std::string data : {
    "dimer_harmonic",
    "trimer_rigid_angle",
    "trimer_harmonic",
    "tetramer_harmonic_no_dihedral",
    "tetramer_harmonic_rigid_bond_angle",
    "tetramer_rigid",
    "tetramer_harmonic",
    "tetramer_branched",
    "pentamer_harmonic",
    }) {
    //for (const std::string num_steps : {"1"}) {
    for (const std::string num_steps : {"1", "4"}) {
      std::string ref = "-1";
      if (num_steps == "4") ref = "0";
      INFO("data " << data << " num_steps " << num_steps << " ref " << ref);
      MonteCarlo mc;
      //mc.set(MakeRandomMT19937({{"seed", "123"}}));
      mc.add(MakeConfiguration({
        {"particle_type0", "../plugin/chain/test/data/" + data + ".fstprt"},
        {"add_particles_of_type0", "1"},
        {"cubic_box_length", "10"}}));
      mc.add(MakePotential(MakeDontVisitModel()));
      mc.add_to_reference(MakePotential(MakeDontVisitModel()));
      mc.set(MakeThermoParams({{"beta", "1"}}));
      mc.set(MakeMetropolis());
      DEBUG("initial energy " << mc.criteria().current_energy());
      if (data == "dimer_harmonic") {
        mc.add(MakeTrialGrow({
          {{"particle_type", "0"}, {"bond", "1"}, {"mobile_site", "1"}, {"anchor_site", "0"}},
        }, {{"num_steps", num_steps}, {"reference_index", ref}}));
      } else if (data == "trimer_rigid_angle" || data == "trimer_harmonic") {
        mc.add(MakeTrialGrow({
          {{"particle_type", "0"}, {"bond", "1"}, {"mobile_site", "1"}, {"anchor_site", "0"}},
          {{"angle", "1"}, {"mobile_site", "2"}, {"anchor_site", "1"}, {"anchor_site2", "0"}},
        }, {{"num_steps", num_steps}, {"reference_index", ref}}));
      } else if (data == "tetramer_harmonic_no_dihedral") {
        mc.add(MakeTrialGrow({
          {{"particle_type", "0"}, {"bond", "1"}, {"mobile_site", "1"}, {"anchor_site", "0"}},
          {{"angle", "1"}, {"mobile_site", "2"}, {"anchor_site", "1"}, {"anchor_site2", "0"}},
          {{"angle", "1"}, {"mobile_site", "3"}, {"anchor_site", "2"}, {"anchor_site2", "1"}},
        }, {{"num_steps", num_steps}, {"reference_index", ref}}));
      } else if (data == "tetramer_harmonic" ||
                 data == "tetramer_harmonic_rigid_bond_angle" ||
                 data == "tetramer_rigid") {
        mc.add(MakeTrialGrow({
          {{"particle_type", "0"}, {"bond", "1"}, {"mobile_site", "1"}, {"anchor_site", "0"}},
          {{"angle", "1"}, {"mobile_site", "2"}, {"anchor_site", "1"}, {"anchor_site2", "0"}},
          {{"dihedral", "1"}, {"mobile_site", "3"}, {"anchor_site", "2"}, {"anchor_site2", "1"}, {"anchor_site3", "0"}},
        }, {{"num_steps", num_steps}, {"reference_index", ref}}));
      } else if (data == "tetramer_branched") {
        mc.add(MakeTrialGrow({
          {{"particle_type", "0"}, {"bond", "1"}, {"mobile_site", "1"}, {"anchor_site", "0"}},
          {{"branch", "1"}, {"mobile_site", "2"}, {"mobile_site2", "3"}, {"anchor_site", "0"}, {"anchor_site2", "1"}},
        }, {{"num_steps", num_steps}, {"reference_index", ref}}));
      } else if (data == "pentamer_harmonic") {
        mc.add(MakeTrialGrow({
          {{"particle_type", "0"}, {"bond", "1"}, {"mobile_site", "1"}, {"anchor_site", "0"}},
          {{"angle", "1"}, {"mobile_site", "2"}, {"anchor_site", "1"}, {"anchor_site2", "0"}},
          {{"dihedral", "1"}, {"mobile_site", "3"}, {"anchor_site", "2"}, {"anchor_site2", "1"}, {"anchor_site3", "0"}},
          {{"dihedral", "1"}, {"mobile_site", "4"}, {"anchor_site", "3"}, {"anchor_site2", "2"}, {"anchor_site3", "1"}},
        }, {{"num_steps", num_steps}, {"reference_index", ref}}));
      }
      //const std::string steps_per = "1";
      const std::string steps_per = "1e3";
      mc.add(MakeLogAndMovie({{"steps_per", steps_per}, {"file_name", "tmp/harmonic"}}));
      mc.add(MakeCheckEnergy({{"steps_per", steps_per}}));
      auto en = MakeEnergy({{"steps_per_write", steps_per}});
      mc.add(en);
      auto bonds = MakeAnalyzeBonds({{"bond_bin_width", "0.05"}});
      // auto bonds = MakeAnalyzeBonds({{"bond_bin_width", "0.05"}, {"steps_per", "1e3"}});
      mc.add(bonds);
      mc.attempt(1e6);
      //INFO(bonds->bond_hist(0).str());
      //INFO(bonds->bond(0).average() << " +/- " << 3*bonds->bond(0).block_stdev());
      DEBUG(bonds->bond(0).str());
      double z_fac = 30;
      if (num_steps != "1") {
        z_fac *= 3;
        if (data == "tetramer_branched") { z_fac *= 2; }
      }
      double l_expect = 1.00167;
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
}

TEST(MonteCarlo, single_butane) {
  MonteCarlo mc;
  //mc.set(MakeRandomMT19937({{"seed", "123"}}));
  mc.add(MakeConfiguration({
    {"particle_type0", "../forcefield/n-butane.fstprt"},
    {"add_particles_of_type0", "1"},
    {"cubic_box_length", "100"}}));
  mc.add(MakePotential(MakeLennardJones(),
                       MakeVisitModelIntra({{"cutoff",  "4"}})));
  mc.set(MakeThermoParams({{"beta", "1"}}));
  mc.set(MakeMetropolis());
  DEBUG("initial energy " << mc.criteria().current_energy());
  mc.add(MakeTrialGrow({
    {{"particle_type", "0"}, {"bond", "1"}, {"mobile_site", "1"}, {"anchor_site", "0"}},
    {{"angle", "1"}, {"mobile_site", "2"}, {"anchor_site", "1"}, {"anchor_site2", "0"}},
    {{"dihedral", "1"}, {"mobile_site", "3"}, {"anchor_site", "2"}, {"anchor_site2", "1"}, {"anchor_site3", "0"}}}));
  //const std::string steps_per = "1";
  const std::string steps_per = "1e3";
  mc.add(MakeLogAndMovie({{"steps_per", steps_per}, {"file_name", "tmp/butane"}}));
  mc.add(MakeCheckEnergy({{"steps_per", steps_per}}));
  auto en = MakeEnergy({{"steps_per_write", steps_per}});
  mc.add(en);
  const int bins = 20;
  auto bonds = MakeAnalyzeBonds({
      {"bond_bin_width", "0.05"}, {"angle_bin_width", "0.05"},
      {"dihedral_bin_width", str(PI/bins)}, {"dihedral_bin_center", str(PI/bins/2)}});
  // auto bonds = MakeAnalyzeBonds({{"bond_bin_width", "0.05"}, {"steps_per", "1e3"}});
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

}  // namespace feasst
