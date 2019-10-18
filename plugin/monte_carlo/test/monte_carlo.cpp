#include <memory>
#include "utils/test/utils.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/test/monte_carlo_test.h"
#include "monte_carlo/include/criteria_metropolis.h"
#include "utils/include/utils_io.h"
#include "math/include/accumulator.h"
#include "system/test/system_test.h"
#include "system/include/long_range_corrections.h"
#include "system/include/visit_model_intra.h"
#include "system/include/visit_model_cell.h"
#include "steppers/include/num_particles.h"

namespace feasst {

TEST(MonteCarlo, serialize) {
  MonteCarlo mc;
  mc_lj(&mc);
//  { std::stringstream ss;
//    mc.serialize(ss);
//    INFO(ss.str());
//    MonteCarlo mc2(ss);
//  }
  MonteCarlo mc2 = test_serialize(mc);
}

TEST(MonteCarlo, NVT_benchmark) {
  MonteCarlo mc;
  mc_lj(&mc);
  mc.set(MakeRandomMT19937({{"seed", "default"}}));
  mc.seek_num_particles(50);
  // mc.seek_num_particles(250);
  // mc.attempt(1e6);  // ~3.5 seconds (now 4.1) with 50
  // mc.seek_num_particles(450);
  // mc.attempt(1e5);  // 15 sec with 450 on slow computer
  mc.attempt(1e3);
  // DEBUG("\n" << mc.timer_str());
}

TEST(MonteCarlo, NVT_SRSW) {
  MonteCarlo mc;
  mc_lj(&mc);
  const int nMol = 500;
  const double rho = 1e-3, length = pow(static_cast<double>(nMol)/rho, 1./3.);
  mc.get_system()->get_configuration()->set_side_length(
    Position().set_vector({length, length, length}));
  mc.seek_num_particles(nMol);
  Accumulator pe;
  const int num_trials = 1e3;
  for (int trial = 0; trial < num_trials; ++trial) {
    mc.attempt(1);  // ~4 seconds
    pe.accumulate(mc.criteria()->current_energy());
  }
  // HWH temperature not set
  DEBUG("pe " << pe.average());
}

TEST(MonteCarlo, GCMC) {
  MonteCarlo mc;
  mc_lj(&mc);
  mc.set(MakeRandomMT19937({{"seed", "default"}}));
  mc.set(MakeCriteriaMetropolis({{"beta", "1.2"}, {"chemical_potential", "-2"}}));
  add_trial_transfer(&mc, {{"particle_type", "0"}});
  //mc.add(MakeTrialAdd({{"particle_type", "0"}}));
  //mc.add(MakeTrialRemove());
  mc.add(MakeNumParticles({{"steps_per_write", "10000"}}));
  // mc.add(MakeTrialTransfer());
//  mc.add(MakeCheckpoint({{"file_name", "tmp/gcmc"}, {"num_hours", "1"}}));
//  for (int i = 0; i < 1e6; ++i) {
//    if (i%100==0) {
//      INFO(mc.system().configuration().num_particles());
//    }
//    mc.attempt(1);
//  }
  mc.attempt(1e4);  // ~4.7 seconds with ~100 particles
}

// // HWH delete
// TEST(MonteCarlo, new_trial) {
//   MonteCarlo mc = mc_lj();
//   {
//     Potential lj_dual_cut(MakeModelLJ(), MakeVisitModelCell());
//     lj_dual_cut.set_model_params(mc.system().configuration());
//     lj_dual_cut.set_model_param("cutoff", 0, 1);
//     mc.add_to_reference(lj_dual_cut);
//   }
//   mc.seek_num_particles(50);
//   mc.set(MakeCriteriaMetropolis({{"beta", "1.2"}, {"chemical_potential", "-4"}}));
//   auto translate = MakeTrialTranslate({{"reference_index", "-1"}, {"num_steps", "1"}});
//   auto add = MakeTrialAdd();
//   //auto add = MakeTrialAdd({{"reference_index", "0"}, {"num_steps", "3"}});
//   auto remove = MakeTrialRemove();
//   //auto remove = MakeTrialRemove({{"reference_index", "0"}, {"num_steps", "3"}});
//   EXPECT_EQ(translate->stages()[0]->rosenbluth().num(), 1);
//   EXPECT_EQ(translate->stages()[0]->reference(), -1);
//   const double old_energy = mc.get_system()->energy();
//   EXPECT_NEAR(old_energy, mc.criteria()->current_energy(), 50*NEAR_ZERO);
//   for (int step = 0; step < 1e4; ++step) {
//     DEBUG("translating");
//     translate->attempt(mc.get_criteria(), mc.get_system());
//     DEBUG("adding");
//     add->attempt(mc.get_criteria(), mc.get_system());
//     DEBUG("removing");
//     remove->attempt(mc.get_criteria(), mc.get_system());
//     DEBUG("num " << mc.system().configuration().num_particles());
//   }
//   const double new_energy = mc.get_system()->energy();
//   EXPECT_NE(old_energy, new_energy);
//   EXPECT_NEAR(new_energy, mc.criteria()->current_energy(), 1e-10);
//   // EXPECT_NEAR(new_energy, mc.criteria()->current_energy(), 500*NEAR_ZERO);
//   EXPECT_TRUE(translate->num_success() > 0);
//   EXPECT_TRUE(add->num_success() > 0);
//   EXPECT_TRUE(remove->num_success() > 0);
//
//   // test PerturbRemove
//   PerturbRemove perturb_remove;
//   auto sel = std::make_shared<TrialSelectParticleOfType>();
//   sel->select(mc.get_system());
//   const int num = mc.system().configuration().num_particles();
//   perturb_remove.perturb(mc.get_system(), sel.get());
//   EXPECT_EQ(num, mc.system().configuration().num_particles());
//   perturb_remove.finalize(mc.get_system());
//   EXPECT_EQ(num - 1, mc.system().configuration().num_particles());
//
//   //DEBUG(mc.timer().str());
// }

TEST(MonteCarlo, grow) {
  for (int i = 0; i < 1; ++i) { // lj dimer
  // for (int i = 1; i < 2; ++i) { // spce
  //for (int i = 0; i < 2; ++i) { // both
    double box_length = 8.;
    std::string data = "../forcefield/data.dimer";
    if (i == 1) {
      box_length=20;
      data = "../forcefield/data.spce";
    }
    MonteCarlo mc;
    mc_lj(&mc, box_length, data);
    mc.add(build_(0, data));  // 0: move
    mc.add(build_(1, data));  // 1: add
    mc.add(build_(2, data));  // 2: remove
    EXPECT_FALSE(mc.trial(0)->stage(0)->trial_select()->is_ghost());  // translate
    EXPECT_FALSE(mc.trial(1)->stage(0)->trial_select()->is_ghost());  // grow
    EXPECT_TRUE (mc.trial(2)->stage(0)->trial_select()->is_ghost());
    EXPECT_FALSE(mc.trial(3)->stage(0)->trial_select()->is_ghost());
    mc.seek_num_particles(3);
    mc.set(MakeCriteriaMetropolis({{"beta", "1.2"}, {"chemical_potential", "-700"}}));
    mc.add(MakeMovie(
     {{"steps_per", "1"},
      {"file_name", "tmp/grow.xyz"}}));
    for (int i = 0; i < 2e1; ++i) {
      mc.attempt(1);
      //mc.configuration().check();
    }
    EXPECT_LT(mc.configuration().num_particles(), 3);
    mc.set(MakeCriteriaMetropolis({{"beta", "1.2"}, {"chemical_potential", "100"}}));
    mc.attempt(2e1);
    EXPECT_GE(mc.configuration().num_particles(), 1);
    mc.configuration().check();
    // INFO(mc.trial(1)->accept().perturbed().str());
  }
}

}  // namespace feasst
