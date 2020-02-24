#include "utils/test/utils.h"
#include "system/include/lennard_jones.h"
#include "math/include/random_mt19937.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/monte_carlo.h"
#include "chain/include/trial_deprotonation.h"
#include "chain/include/trial_protonation.h"
#include "chain/include/trial_swap_sites.h"
#include "chain/include/utils_chain.h"

namespace feasst {

TEST(TrialDeprotonation, deprotonationANDprotonation) {
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
//  random = MakeRandomMT19937({{"seed", "1581704232"}});
  auto criteria = MakeMetropolis({
    {"beta", "1.0"},
    {"chemical_potential0", "1"},
    {"chemical_potential1", "-5"},
    {"pH", "1"},
  });
  auto trial = MakeTrialDeprotonation({
    {"reactant_type", "0"},
    {"reactant_site_type", "0"},
    {"new_site_type", "1"},
    {"add_type", "1"},
  });
  trial->precompute(criteria.get(), &system);
  while (trial->num_success() == 0) {
    trial->attempt(criteria.get(), &system, random.get());
    ASSERT(trial->num_attempts() < 1e3, "not accepting");
  }
  EXPECT_EQ(trial->num_success(), 1);
  EXPECT_EQ(system.configuration().num_particles(), 2);
  const int site_index = trial->stage(0)->trial_select()->mobile().site_index(0, 0);
  EXPECT_TRUE(site_index == 1 ||
              site_index == 3 ||
              site_index == 4 ||
              site_index == 5 ||
              site_index == 6 ||
              site_index == 7 ||
              site_index == 8 ||
              site_index == 9);
  EXPECT_EQ(system.configuration().particle(0).site(site_index).type(), 1);

  // attempt a protonation, now that it has been deprotonated
  auto prot = MakeTrialProtonation({
    {"reactant_type", "0"},
    {"reactant_site_type", "1"},
    {"new_site_type", "0"},
    {"remove_type", "1"},
  });
  prot->precompute(criteria.get(), &system);
  while (prot->num_success() == 0) {
    prot->attempt(criteria.get(), &system, random.get());
    ASSERT(prot->num_attempts() < 1e3, "not accepting");
  }
  EXPECT_EQ(prot->num_success(), 1);
  EXPECT_EQ(system.configuration().num_particles(), 1);
  EXPECT_EQ(system.configuration().particle(0).type(), 0);
  const int site_index2 = prot->stage(0)->trial_select()->mobile().site_index(0, 0);
  EXPECT_TRUE(site_index2 == 0 || site_index2 == site_index);
  EXPECT_EQ(system.configuration().particle(0).site(site_index2).type(), 0);

  // attempt to swap sites
  auto swap = MakeTrialSwapSites({
    {"particle_type", "0"},
    {"site_type1", "0"},
    {"site_type2", "1"},
  });
  swap->precompute(criteria.get(), &system);
  while (swap->num_success() == 0) {
    swap->attempt(criteria.get(), &system, random.get());
    ASSERT(swap->num_attempts() < 1e3, "not accepting");
  }
  EXPECT_EQ(swap->num_success(), 1);
  const int site_index3 = swap->stage(0)->trial_select()->mobile().site_index(0, 0);
  const int site_index4 = swap->stage(1)->trial_select()->mobile().site_index(0, 0);
  EXPECT_EQ(1, system.configuration().particle(0).site(site_index3).type());
  EXPECT_EQ(0, system.configuration().particle(0).site(site_index4).type());

  MonteCarlo monte_carlo;
  monte_carlo.set(system);
  monte_carlo.add(criteria);
  add_deprotonation(&monte_carlo, {
    {"reactant_type", "0"},
    {"reactant_site_type", "0"},
    {"new_site_type", "1"},
    {"add_type", "1"},
  });
  while (monte_carlo.trial(0)->num_success() == 0 ||
         monte_carlo.trial(1)->num_success() == 0) {
    monte_carlo.attempt();
    ASSERT(monte_carlo.trials().num_attempts() < 1e4, "not accepting");
  }
}

}  // namespace feasst
