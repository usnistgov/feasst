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

namespace feasst {

TEST(MonteCarlo, serialize) {
  MonteCarlo mc = mc_lj(), mc2 = test_serialize(mc);
}

TEST(MonteCarlo, NVT_benchmark) {
  seed_random_by_date();
  seed_random();
  MonteCarlo mc = mc_lj();
  mc.seek_num_particles(50);
  // mc.attempt(1e6);  // ~3.5 seconds
  mc.attempt(1e4);
  // DEBUG("\n" << mc.timer_str());
}

TEST(MonteCarlo, NVT_SRSW) {
  seed_random_by_date();
  MonteCarlo mc = mc_lj();
  const int nMol = 500;
  const double rho = 1e-3, length = pow(static_cast<double>(nMol)/rho, 1./3.);
  mc.get_system()->get_configuration()->set_side_length(
    Position().set_vector({length, length, length}));
  mc.seek_num_particles(nMol);
  Accumulator pe;
  for (int trial = 0; trial < 1e3; ++trial) {
    mc.attempt(1);  // ~4 seconds
    pe.accumulate(mc.criteria()->current_energy());
  }
  // HWH temperature not set
  DEBUG("pe " << pe.average());
}

TEST(MonteCarlo, GCMC) {
  MonteCarlo mc = mc_lj();
  mc.add(MakeTrialAdd());
  mc.add(MakeTrialRemove());
  // mc.add(MakeTrialTransfer());
  mc.add(MakeCheckpoint({{"file_name", "tmp/gcmc"}, {"num_hours", "1"}}));
  mc.attempt(1e4);
}

// HWH delete
TEST(MonteCarlo, new_trial) {
  seed_random_by_date();
//  seed_random(1558629935);
  //seed_random(1558643002);
//  seed_random(1558662620);
  seed_random();
  MonteCarlo mc = mc_lj();
  {
    Potential lj_dual_cut(MakeModelLJ(), MakeVisitModelCell());
    lj_dual_cut.set_model_params(mc.system().configuration());
    lj_dual_cut.set_model_param("cutoff", 0, 1);
    mc.add_to_reference(lj_dual_cut);
  }
  mc.seek_num_particles(50);
  mc.set(MakeCriteriaMetropolis({{"beta", "1.2"}, {"chemical_potential", "-4"}}));
  auto translate = MakeTrialTranslate({{"reference_index", "-1"}, {"num_steps", "1"}});
  auto add = MakeTrialAdd();
  //auto add = MakeTrialAdd({{"reference_index", "0"}, {"num_steps", "3"}});
  auto remove = MakeTrialRemove();
  //auto remove = MakeTrialRemove({{"reference_index", "0"}, {"num_steps", "3"}});
  EXPECT_EQ(translate->stages()[0]->rosenbluth().num(), 1);
  EXPECT_EQ(translate->stages()[0]->reference(), -1);
  const double old_energy = mc.get_system()->energy();
  EXPECT_NEAR(old_energy, mc.criteria()->current_energy(), 50*NEAR_ZERO);
  for (int step = 0; step < 1e4; ++step) {
    DEBUG("translating");
    translate->attempt(mc.get_criteria(), mc.get_system());
    DEBUG("adding");
    add->attempt(mc.get_criteria(), mc.get_system());
    DEBUG("removing");
    remove->attempt(mc.get_criteria(), mc.get_system());
    DEBUG("num " << mc.system().configuration().num_particles());
  }
  const double new_energy = mc.get_system()->energy();
  EXPECT_NE(old_energy, new_energy);
  EXPECT_NEAR(new_energy, mc.criteria()->current_energy(), 1e-10);
  // EXPECT_NEAR(new_energy, mc.criteria()->current_energy(), 500*NEAR_ZERO);
  EXPECT_TRUE(translate->num_success() > 0);
  EXPECT_TRUE(add->num_success() > 0);
  EXPECT_TRUE(remove->num_success() > 0);

  // test PerturbRemove
  PerturbRemove perturb_remove;
  auto sel = std::make_shared<TrialSelectParticleOfType>();
  sel->select(mc.get_system());
  const int num = mc.system().configuration().num_particles();
  perturb_remove.perturb(mc.get_system(), sel.get());
  EXPECT_EQ(num, mc.system().configuration().num_particles());
  perturb_remove.finalize(mc.get_system());
  EXPECT_EQ(num - 1, mc.system().configuration().num_particles());

  //DEBUG(mc.timer().str());
}

}  // namespace feasst
