#include <memory>
#include <cmath>
#include "utils/test/utils.h"
#include "utils/include/utils_io.h"
#include "utils/include/checkpoint.h"
#include "math/include/accumulator.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/domain.h"
#include "system/include/long_range_corrections.h"
#include "system/include/visit_model_intra.h"
#include "system/include/visit_model_cell.h"
#include "system/include/utils.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/test/monte_carlo_test.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/constrain_num_particles.h"
#include "monte_carlo/include/utils.h"
#include "monte_carlo/include/seek_num_particles.h"
#include "steppers/include/num_particles.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/movie.h"
#include "steppers/include/utils.h"

namespace feasst {

//TEST(MonteCarlo, tutorial) {
//  MonteCarlo mc(
//    MakeRandomMT19937({{"seed", "time"}}),
//    Configuration(feasst::MakeDomain({{"cubic_box_length", "8"}}),
//      {{"particle_type", install_dir() + "/forcefield/data.lj")}}));
//  mc.add(feasst::Potential(feasst::MakeLennardJones()));
//  mc.add(feasst::Potential(feasst::MakeLongRangeCorrections()));
//  mc.add(feasst::MakeMetropolis(
//    {{"beta", "1.2"}, {"chemical_potential", "1."}}));
//  mc.add(feasst::MakeTrialTranslate(
//    {{"tunable_param", "2."}, {"tunable_target_acceptance", "0.2"}}));
//  const int steps_per = 1e3;
//  mc.add(feasst::MakeTuner({{"steps_per", feasst::str(steps_per)}}));
//  mc.seek_num_particles(50);
//  mc.add(feasst::MakeLog({{"steps_per", feasst::str(steps_per)}}));
//  mc.add(feasst::MakeMovie(
//   {{"steps_per", feasst::str(steps_per)}, {"file_name", "movie.xyz"}}));
//  mc.add(feasst::MakeCheckEnergy(
//   {{"steps_per", feasst::str(steps_per)}, {"tolerance", "1e-8"}}));
//  mc.attempt(1e5);
//}

TEST(MonteCarlo, serialize) {
  MonteCarlo mc;
  mc.set(lennard_jones());
  mc.set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc.set(MakeCheckpoint({{"num_hours", "0.0001"}, {"file_name", "tmp/ljrst"}}));
  add_common_steppers(&mc, {{"steps_per", str(1e4)},
                            {"file_append", "tmp/lj"}});
//  { std::stringstream ss;
//    mc.serialize(ss);
//    INFO(ss.str());
//    MonteCarlo mc2(ss);
//  }
  MonteCarlo mc2 = test_serialize(mc);
}

TEST(MonteCarlo, NVT_benchmark) {
  MonteCarlo mc;
  mc.set(lennard_jones());
  mc.set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  add_common_steppers(&mc, {{"steps_per", str(1e4)},
                            {"file_append", "tmp/lj"}});
  //mc.set(0, Potential(MakeLennardJones(),
  //  MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
//  mc.add_to_optimized(Potential(MakeLennardJones(), //HWH: prevents ModelEmpty... how to remove?
//                                MakeVisitModelOptLJ()));
  mc.set(MakeRandomMT19937({{"seed", "default"}}));
  SeekNumParticles(50).with_trial_add().run(&mc);
  // mc.seek_num_particles(250);
  // mc.attempt(1e6);  // 5.4s with 50 (see opt_lj for 4.3s)
  // mc.seek_num_particles(450);
  // mc.attempt(1e5);  // 15 sec with 450 on slow computer
  mc.attempt(1e3);
  // DEBUG("\n" << mc.timer_str());
}

TEST(MonteCarlo, NVT_cell_benchmark) {
  MonteCarlo mc;
  mc.add(Configuration(MakeDomain({{"cubic_box_length", "8"},
                                   {"init_cells", "1"}}),
                       {{"particle_type", "../forcefield/data.lj"}}));
  mc.add(Potential(MakeLennardJones()));
  mc.add_to_reference(Potential(MakeLennardJones(), MakeVisitModelCell()));
  mc.set(MakeRandomMT19937({{"seed", "default"}}));
  mc.set(MakeMetropolis({{"beta", "1.2"}}));
  mc.add(MakeTrialTranslate({
    {"weight", "1"},
    {"reference_index", "0"},
    {"num_steps", "10"},
    {"tunable_param", "1"}}));
//  const int steps_per = 1e3;
//  mc.add(MakeCheckEnergy({{"steps_per", str(steps_per)}, {"tolerance", "1e-10"}}));
  SeekNumParticles(50)
    .with_metropolis({{"beta", "0.1"}, {"chemical_potential", "10"}})
    .with_trial_add()
    .run(&mc);
  add_common_steppers(&mc, {{"steps_per", str(1e4)},
                            {"file_append", "tmp/cell"}});
  //mc.attempt(1e5);  // 3 sec with 50 particles, 10 steps on hwhdesk after site cells integer
  //mc.attempt(1e5);  // 4.4 sec with 50 particles, 10 steps on hwhdesk after cell group opt
  //mc.attempt(1e5);  // 4.5 sec with 50 particles, 10 steps on hwhdesk after domain shift opt
  //mc.attempt(1e5);  // 5.6 sec with 50 particles, 10 steps on hwhdesk
  //mc.attempt(1e5);  // 5.1 sec with 50 particles, 10 steps on hwhdesk after opt cell_id
  mc.attempt(1e3);
  INFO(mc.criteria().current_energy());
  // mc.attempt(1e6);
}

TEST(MonteCarlo, NVT_SRSW) {
  MonteCarlo mc;
  mc.set(lennard_jones());
  mc.set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  add_common_steppers(&mc, {{"steps_per", str(1e4)},
                            {"file_append", "tmp/lj"}});
  const int nMol = 500;
  const double rho = 1e-3, length = std::pow(static_cast<double>(nMol)/rho, 1./3.);
  mc.get_system()->get_configuration()->set_side_lengths(
    Position().set_vector({length, length, length}));
  SeekNumParticles(nMol).with_trial_add().run(&mc);
  Accumulator pe;
  const int num_trials = 1e3;
  for (int trial = 0; trial < num_trials; ++trial) {
    mc.attempt(1);  // ~4 seconds
    pe.accumulate(mc.criteria().current_energy());
  }
  // HWH temperature not set
  DEBUG("pe " << pe.average());
}

TEST(MonteCarlo, GCMC) {
  MonteCarlo mc;
  mc.set(lennard_jones());
  mc.set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  add_common_steppers(&mc, {{"steps_per", str(1e4)},
                            {"file_append", "tmp/lj"}});
  mc.set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "-6"}}));
  add_trial_transfer(&mc, {{"particle_type", "0"}});
  mc.add(MakeNumParticles({{"steps_per_write", str(1e5)},
                           {"file_name", "tmp/ljnum.txt"}}));
  mc.attempt(1e4);
}

TEST(MonteCarlo, GCMC_cell) {
  MonteCarlo mc;
  mc.set(lennard_jones({{"dual_cut", "1."}}));
  mc.set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc.add(MakeTrialTranslate({{"weight", "1."},
                             {"tunable_param", "1."},
                             {"reference_index", "0"},
                             {"num_steps", "4"}}));
  add_common_steppers(&mc, {{"steps_per", str(1e4)},
                            {"file_append", "tmp/lj"}});
  mc.set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "-6"}}));
  add_trial_transfer(&mc,
    { {"particle_type", "0"},
      {"num_steps", "4"},
      {"reference_index", "0"}});
  mc.add(MakeNumParticles({{"steps_per_write", str(1e5)},
                           {"file_name", "tmp/ljnum.txt"}}));
  mc.attempt(1e4);
}

// HWH this test is known to fail infrequently
TEST(MonteCarlo, grow) {
  for (int i = 0; i < 1; ++i) { // lj dimer
  // for (int i = 1; i < 2; ++i) { // spce
  //for (int i = 0; i < 2; ++i) { // both
    double box_length = 8.;
    std::string data = "forcefield/data.dimer";
    if (i == 1) {
      box_length=20;
      data = "forcefield/data.spce";
    }
    MonteCarlo mc;
    mc.set(lennard_jones({{"cubic_box_length", str(box_length)},
                          {"particle", data}}));
    mc.set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "1."}}));
    mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
    add_common_steppers(&mc, {{"steps_per", str(1e4)},
                              {"file_append", "tmp/lj"}});
    mc.add(build_(0, data));  // 0: move
    mc.add(build_(1, data));  // 1: add
    mc.add(build_(2, data));  // 2: remove
    EXPECT_FALSE(mc.trial(0).stage(0).trial_select().is_ghost());  // translate
    EXPECT_FALSE(mc.trial(1).stage(0).trial_select().is_ghost());  // grow
    EXPECT_TRUE (mc.trial(2).stage(0).trial_select().is_ghost());
    EXPECT_FALSE(mc.trial(3).stage(0).trial_select().is_ghost());
    SeekNumParticles(3).with_trial_add().run(&mc);
    mc.set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "-700"}}));
    mc.add(MakeMovie(
     {{"steps_per", "1"},
      {"file_name", "tmp/grow.xyz"}}));
    for (int i = 0; i < 2e1; ++i) {
      mc.attempt(1);
      //mc.configuration().check();
    }
    EXPECT_LT(mc.configuration().num_particles(), 3);
    mc.set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "100"}}));
    mc.attempt(2e1);
    EXPECT_GE(mc.configuration().num_particles(), 1);
    mc.configuration().check();
    // INFO(mc.trial(1)->accept().perturbed().str());
  }
}

TEST(MonteCarlo, ConstrainNumParticles) {
  for (const double minimum : {0, 1}) {
    MonteCarlo monte_carlo;
    monte_carlo.set(lennard_jones());
    monte_carlo.set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "1."}}));
    monte_carlo.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
    add_common_steppers(&monte_carlo, {{"steps_per", str(1e4)},
                                       {"file_append", "tmp/lj"}});
    SeekNumParticles(1).with_trial_add().run(&monte_carlo);
    monte_carlo.add(MakeMetropolis(
      MakeConstrainNumParticles({{"minimum", str(minimum)},
                                 {"maximum", str(minimum+1)}}),
      {{"beta", "0.2"}, {"chemical_potential", "-20."}}));
    add_trial_transfer(&monte_carlo, {{"particle_type", "0"}});
    const int index = monte_carlo.num_analyzers();
    monte_carlo.add(MakeNumParticles({{"steps_per_write", "10000"}}));
    monte_carlo.attempt(14);
    if (minimum == 0) {
      EXPECT_LE(monte_carlo.analyze(index).accumulator().average(), 1);
    } else if (minimum == 1) {
      EXPECT_GE(monte_carlo.analyze(index).accumulator().average(), 1);
    }
  }
}

}  // namespace feasst
