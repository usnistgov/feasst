#include "utils/test/utils.h"
#include "system/include/lennard_jones.h"
#include "math/include/random_mt19937.h"
#include "monte_carlo/include/metropolis.h"
#include "chain/include/trial_synthesis.h"

namespace feasst {

TEST(TrialSynthesis, synthesis) {
  System system;
  {
    Configuration config({
      {"particle_type0", "../plugin/chain/forcefield/data.chain10titratable"},
      {"particle_type1", "../forcefield/data.lj"},
      {"cubic_box_length", "20"}});
    config.add_particle_of_type(0);
    system.add(config);
  }
  system.add(Potential(MakeLennardJones()));
  system.energy();
  auto random = MakeRandomMT19937();
  auto criteria = MakeMetropolis({
    {"beta", "1.0"},
    {"chemical_potential0", "1"},
    {"chemical_potential1", "1"},
  });
  auto trial = MakeTrialSynthesis({
    {"reactant_type", "0"},
    {"reactant_site_type", "0"},
    {"new_site_type", "1"},
    {"product_type", "1"},
  });
  trial->precompute(criteria.get(), &system);
  while (trial->num_success() == 0) {
    trial->attempt(criteria.get(), &system, random.get());
    ASSERT(trial->num_attempts() < 1e3, "not accepting");
  }
  EXPECT_EQ(trial->num_success(), 1);
  EXPECT_EQ(system.configuration().num_particles(), 2);
  const int site_index = trial->stage(0)->trial_select()->mobile().site_index(0, 0);
  EXPECT_TRUE(site_index == 0 ||
              site_index == 3 ||
              site_index == 4 ||
              site_index == 5 ||
              site_index == 6 ||
              site_index == 7 ||
              site_index == 8 ||
              site_index == 9);
  EXPECT_EQ(system.configuration().particle(0).site(site_index).type(), 1);
}

}  // namespace feasst
