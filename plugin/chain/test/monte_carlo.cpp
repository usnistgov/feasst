#include <memory>
#include "utils/test/utils.h"
#include "utils/include/io.h"
#include "math/include/accumulator.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/domain.h"
#include "system/include/ideal_gas.h"
#include "system/include/hard_sphere.h"
#include "system/include/utils.h"
#include "system/include/visit_model_intra_map.h"
#include "models/include/square_well.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trials.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/seek_num_particles.h"
#include "monte_carlo/include/trial_compute_move.h"
#include "steppers/include/log.h"
#include "steppers/include/check_energy_and_tune.h"
#include "steppers/include/log_and_movie.h"
#include "steppers/include/movie.h"
#include "steppers/include/tuner.h"
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
#include "chain/include/check_rigid_bonds.h"
#include "chain/include/trials.h"
#include "chain/include/trial_grow.h"
#include "chain/include/trial_grow_linear.h"
#include "chain/include/analyze_bonds.h"
#include "chain/test/system_chain.h"

namespace feasst {

TEST(MonteCarlo, chain) {
  MonteCarlo mc;
  //mc.set(MakeRandomMT19937({{"seed", "1602167050"}}));
  mc.set(chain_system());
  mc.set(MakeThermoParams({{"beta", "1"}, {"chemical_potential", "1."}}));
  mc.set(MakeMetropolis());
  SeekNumParticles(1).with_trial_add().run(&mc);
  mc.add(MakeTrialTranslate({
    {"weight", "1."},
    {"tunable_param", "1."},
    {"reference_index", "0"},
    {"num_steps", "2"},
  }));
  mc.add(MakeTrialRotate({
    {"weight", "1."},
    {"tunable_param", "20."},
    {"reference_index", "0"},
    {"num_steps", "2"},
  }));
  mc.add(MakeTrialPivot({
    {"weight", "1."},
    {"tunable_param", "20."},
    {"max_length", "30"},
    {"reference_index", "0"},
    {"num_steps", "2"},
  }));
  mc.add(MakeTrialReptate({
    {"weight", "1."},
    {"max_length", "1"},
//    {"reference_index", "0"},
//    {"num_steps", "2"},
  }));
  mc.add(MakeTrialCrankshaft({
    {"weight", "1."},
    {"tunable_param", "25."},
    {"max_length", "5."},
    {"reference_index", "0"},
    {"num_steps", "2"},
  }));
  mc.add(MakeTrialGrowLinear(
    MakeTrialComputeMove(),
    {
//      {"weight", "0.1"},
      {"particle_type", "0"},
      {"num_steps", "3"},
      {"reference_index", "0"},
    }
  ));
  const int steps_per = 1e2;
  mc.add(MakeLog({
    {"steps_per", str(steps_per)},
    {"file_name", "tmp/chainlog.txt"},
  }));
  mc.add(MakeMovie({
    {"steps_per", str(steps_per)},
    {"file_name", "tmp/chain10movie.xyz"},
  }));
  mc.add(MakeCheckEnergy({
    {"steps_per", str(steps_per)},
    {"tolerance", "1e-10"},
  }));
  mc.add(MakeTuner({{"steps_per", str(steps_per)}}));
  mc.add(MakeCheckRigidBonds({{"steps_per", str(steps_per)}}));
  mc.attempt(3e2);

  MonteCarlo mc2 = test_serialize(mc);
  EXPECT_EQ(mc2.analyzers().size(), 3);
}

// HWH this test is known to fail infrequently
TEST(MonteCarlo, TrialGrow) {
  for (const std::string particle : {"lj", "spce"}) {
    double box_length = 8.;
    std::string data = "forcefield/data.dimer";
    if (particle == "spce") {
      box_length=20;
      data = "forcefield/data.spce";
    }
    MonteCarlo mc;
    mc.set(MakeRandomMT19937({{"seed", "123"}}));
    mc.set(lennard_jones({{"cubic_box_length", str(box_length)},
                          {"particle", data}}));
    mc.set(MakeThermoParams({{"beta", "1.2"}, {"chemical_potential", "-700"}}));
    mc.set(MakeMetropolis());
    SeekNumParticles(3)
      .with_thermo_params({{"beta", "1.2"}, {"chemical_potential", "1."}})
      .with_trial_add().run(&mc);
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
    mc.add(MakeTrialGrow(grow_args, {{"num_steps", "4"}}));
    EXPECT_EQ(4, mc.trial(0).stage(0).rosenbluth().num());
    EXPECT_EQ(5, mc.trial(0).stage(1).rosenbluth().num());
    mc.add(MakeLogAndMovie({{"steps_per", str(1e0)}, {"file_name", "tmp/lj"}}));
    mc.add(MakeCheckEnergyAndTune({{"steps_per", str(1e0)}, {"tolerance", str(1e-9)}}));
    mc.add(MakeCheckRigidBonds({{"steps_per", str(1e0)}}));
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
  MonteCarlo mc;
  //mc.set(MakeRandomMT19937({{"seed", "123"}}));
  mc.add(Configuration(MakeDomain({{"cubic_box_length", "30"}}),
    {{"particle_type", "../plugin/chain/forcefield/" + data},
     {"set_cutoff_min_to_sigma", "true"}}));
  mc.add(MakePotential(MakeHardSphere()));
  if (is_found_in(data, "fullangflex")) {
    mc.add(MakePotential(MakeHardSphere(),
                     MakeVisitModelIntraMap({{"exclude_bonds", "true"}})));
  }
  mc.set(MakeThermoParams({{"beta", "1."}, {"chemical_potential", "1"}}));
  mc.set(MakeMetropolis());
  SeekNumParticles(num).with_trial_add().run(&mc);
  //SeekNumParticles(10).with_trial_add().run(&mc);
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
  if (!is_found_in(data, "flex")) {
    mc.add(MakeCheckRigidBonds({{"steps_per", str(steps_per)}}));
  }
  return mc;
}

TEST(MonteCarlo, cg7mab2) {
  for (const std::string data : {
      "data.cg7mab2",
      "data.cg7mab2closed",
      "data.cg7mab2extended",
      "data.cg7mab2flex",
      "data.cg7mab2fullangflex",
    }) {
    cg7mab2(data, 1).attempt(1e2);
  }
}

TEST(MonteCarlo, cg7mab2_LONG) {
  for (const std::string data : {
      "data.cg7mab2",
      "data.cg7mab2closed",
      "data.cg7mab2extended",
      "data.cg7mab2flex",
      "data.cg7mab2fullangflex",
    }) {
    cg7mab2(data, 10, 1e4).attempt(1e6);
  }
}

MonteCarlo test_avb(const bool avb2, const bool avb4 = true) {
  MonteCarlo mc;
  mc.set(MakeRandomMT19937({{"seed", "time"}}));
  //mc.set(MakeRandomMT19937({{"seed", "123"}}));
  Configuration config(MakeDomain({{"side_length0", "6"}, {"side_length1", "6"}}),
    {{"particle_type", "../plugin/chain/forcefield/data.heterotrimer2d"}});
  mc.add(config);
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
  SeekNumParticles(10).with_trial_add().run(&mc);
  const std::string steps_per = feasst::str(1e5);
  mc.add(MakeEnergy());
  mc.add(MakeLogAndMovie({{"file_name", "tmp/trimer2d"}, {"steps_per", steps_per}}));
  mc.add(MakeChirality2D());
  mc.add(MakeAnalyzeBonds());
  mc.add(MakeCheckEnergyAndTune({{"steps_per", steps_per}}));
  MonteCarlo mc2 = test_serialize(mc);
  mc2.attempt(1e7);
  const Analyze& chiral = SeekAnalyze().reference("Chirality2D", mc2);
  EXPECT_NEAR(chiral.accumulator().average(), 10., NEAR_ZERO);
  EXPECT_NEAR(chiral.accumulator().stdev(), 0., NEAR_ZERO);
  return mc2;
}

const double z_factor = 10.;

TEST(MonteCarlo, heterotrimer2d_LONG) {
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
  mc.set(MakeRandomMT19937({{"seed", "123"}}));
  //mc.set(MakeRandomMT19937({{"seed", "1610132694"}}));
  mc.set(lennard_jones({{"particle", "forcefield/data.dimer"},
                        {"cubic_box_length", "6"},
                        {"lrc", "false"}}));
  SeekNumParticles(5)
    .with_thermo_params({{"beta", "1"}, {"chemical_potential", "1"}})
    .with_metropolis()
    .with_trial_add()
    .run(&mc);
  auto neigh = MakeEnergyMapNeighbor();
  mc.set(0, MakePotential(MakeLennardJones(), MakeVisitModel(MakeVisitModelInner(neigh))));
  mc.add_to_reference(MakePotential(MakeLennardJones()));
  mc.set(MakeThermoParams({{"beta", "1"}, {"chemical_potential", "1"}}));
  mc.set(MakeMetropolis());
//  mc.add(MakeTrialTranslate());
//  mc.add(MakeTrialRotate({{"tunable_param", "50"}}));
  mc.add(MakeTrialGrow({
    {{"regrow", "true"}, {"particle_type", "0"}, {"site", "0"}},
    {{"bond", "true"}, {"mobile_site", "1"}, {"anchor_site", "0"}}},
    {{"num_steps", "4"}, {"reference_index", "0"}}));
//  mc.add(MakeTrialGrow({
//    {{"regrow", "true"}, {"particle_type", "0"}, {"site", "1"}},
//    {{"bond", "true"}, {"mobile_site", "0"}, {"anchor_site", "1"}}},
//    {{"num_steps", "4"}}));
  mc.add(MakeLogAndMovie({{"steps_per", "100"}, {"file_name", "tmp/dimer"}}));
  for (int i = 0; i < 1e1; ++i) {
    mc.attempt(1);
    neigh->check();
  }
  EXPECT_NEAR(mc.criteria().current_energy(), neigh->total_energy(), 1e-12);
}

}  // namespace feasst
