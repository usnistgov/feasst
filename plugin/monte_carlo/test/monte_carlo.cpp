#include <memory>
#include "utils/test/utils.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/test/monte_carlo_test.h"
#include "monte_carlo/include/metropolis.h"
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
  //mc.set(0, Potential(MakeLennardJones(),
  //  MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
//  mc.add_to_optimized(Potential(MakeLennardJones(), //HWH: prevents ModelEmpty... how to remove?
//                                MakeVisitModelOptLJ()));
  mc.set(MakeRandomMT19937({{"seed", "default"}}));
  mc.seek_num_particles(50);
  // mc.seek_num_particles(250);
  //mc.attempt(1e6);  // 5.4s with 50 (see opt_lj for 4.3s)
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
  mc.set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "-6"}}));
  add_trial_transfer(&mc, {{"particle_type", "0"}});
  mc.add(MakeNumParticles({{"steps_per_write", "1000"},
                           {"file_name", "tmp/ljnum.txt"}}));
  mc.attempt(1e4);
}

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

}  // namespace feasst
