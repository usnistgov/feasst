#include "utils/test/utils.h"
#include "utils/include/utils.h"
#include "system/include/visit_model.h"
#include "math/include/constants.h"
#include "math/include/random_mt19937.h"
#include "system/test/system_test.h"
#include "cluster/include/energy_map_all.h"

namespace feasst {

TEST(EnergyMap, energy_map) {
  Configuration config = lj_sample();
  LennardJones model;
  VisitModel visit(MakeVisitModelInner(MakeEnergyMapAll()));
  visit.precompute(&config);
  model.compute(&config, &visit);
  // HWH: perhaps figure out how to remove this finalize
  visit.finalize(config.selection_of_all());
  const double en_lj_expect = -16.790321304625856;
  EXPECT_NEAR(en_lj_expect, visit.energy(), NEAR_ZERO);
  EXPECT_NEAR(en_lj_expect,
              visit.inner()->energy_map()->total_energy(),
              1e-13);

  // find neighbors within 3 of first particle manually
  std::vector<int> neighs, neighs_known = {1, 4, 5, 8, 9, 11, 13, 14, 15, 21, 25, 27, 28};
  for (int ipart = 1; ipart < config.num_particles(); ++ipart) {
    Position pos = config.particle(ipart).site(0).position();
    pos.subtract(config.particle(0).site(0).position());
    if (pos.distance() < 3) {
      neighs.push_back(ipart);
    }
  }
  // INFO("num " << neighs.size());
  // INFO(feasst_str(neighs));
  EXPECT_EQ(neighs, neighs_known);

  RandomMT19937 random;
  Select neighs2;
  visit.inner()->energy_map()->neighbors(
    MakeNeighborCriteria().get(),
    config,
    0, 0, 0,
    &random,
    &neighs2);
  const int neighbor = random.const_element(neighs2.particle_indices());
  EXPECT_EQ(13, static_cast<int>(neighs2.num_sites()));
  EXPECT_TRUE(find_in_list(neighbor, neighs2.particle_indices()));
}

}  // namespace feasst
