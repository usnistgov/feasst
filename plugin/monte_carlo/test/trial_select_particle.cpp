#include "utils/test/utils.h"
#include "monte_carlo/test/monte_carlo_utils.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"
#include "monte_carlo/include/trial_select_particle.h"

namespace feasst {

TEST(TrialSelectParticle, serialize) {
  TrialSelectParticle add;
//  std::stringstream ss;
//  add.serialize(ss);
//  DEBUG(ss.str());
//  TrialSelectParticle add2(ss);
//  add2.serialize(ss);
//  DEBUG(ss.str());
  TrialSelectParticle add2 = test_serialize(add);
}

TEST(TrialSelectParticle, exclude_perturbed) {
  System system;
  { auto config = MakeConfiguration({
      {"particle_type0", "../particle/lj.txt"},
      {"particle_type1", "../particle/atom.txt"}});
    config->add_particle_of_type(0);
    config->add_particle_of_type(0);
    config->add_particle_of_type(1);
    system.add(config);
  }
  auto sel1 = MakeTrialSelectParticle({{"particle_type", "0"}});
  auto sel2 = MakeTrialSelectParticle({{"particle_type", "0"},
                                       {"exclude_perturbed", "true"}});
  sel1->precompute(&system);
  sel2->precompute(&system);
  auto ran = MakeRandomMT19937();
  ran = MakeRandomMT19937({{"seed", "1583520072"}});
  sel1->sel(&system, ran.get());
  const int part1_index = sel1->mobile().particle_index(0);
  EXPECT_TRUE(part1_index == 0 || part1_index == 1);
  DEBUG(sel1->mobile().str());
  sel2->select(sel1->mobile(), &system, ran.get(), NULL);
  if (part1_index == 0) {
    EXPECT_EQ(sel2->mobile().particle_index(0), 1);
  } else if (part1_index == 1) {
    EXPECT_EQ(sel2->mobile().particle_index(0), 0);
  }
  DEBUG(sel2->mobile().str());
  test_serialize(*sel2);

  // now try the same with ghosts
  sel1 = MakeTrialSelectParticle({{"particle_type", "0"}, {"ghost", "true"}});
  sel2 = MakeTrialSelectParticle({{"particle_type", "0"}, {"ghost", "true"},
                                  {"exclude_perturbed", "true"}});
  sel1->precompute(&system);
  sel1->sel(&system, ran.get());
  sel2->precompute(&system);
  sel2->select(sel1->mobile(), &system, ran.get(), NULL);
  EXPECT_NE(sel1->mobile().particle_index(0), sel2->mobile().particle_index(0));
}

}  // namespace feasst
