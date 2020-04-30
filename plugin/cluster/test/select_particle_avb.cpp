#include "utils/test/utils.h"
#include "cluster/include/select_particle_avb.h"
#include "math/include/random_mt19937.h"
#include "system/include/lennard_jones.h"
#include "monte_carlo/include/perturb_translate.h"
#include "cluster/include/energy_map_all.h"
#include "cluster/include/energy_map_neigh.h"
#include "cluster/include/perturb_rotate_com.h"

namespace feasst {

TEST(SelectParticleAVB, serialize) {
  System system;
  {
    Configuration config(MakeDomain({{"cubic_box_length", "8"}}),
                         {{"particle_type", "../forcefield/data.lj"}});
    config.add_particle_of_type(0);
    config.add_particle_of_type(0);
    config.add_particle_of_type(0);
    config.update_positions({{0, 0, 0},
                             {-1.25, 0, 0},
                             {2.9, 0, 0}});
    system.add(config);
  }
  system.add(Potential(MakeLennardJones(),
                    MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
  system.energy();
  system.finalize();

  auto ran = MakeRandomMT19937();
  //ran = MakeRandomMT19937({{"seed", "1580154124"}});

  for (const double max_dist : {3.}) {
  //for (const double max_dist : {1.5}) {
  //for (const double max_dist : {1.5, 3.}) {
    SelectParticleAVB sel(
      MakeNeighborCriteria({{"maximum_distance", str(max_dist)}}),
      {{"grand_canonical", "true"},
       {"particle_type", "0"}});
    sel.precompute(&system);
    auto sel2 = test_serialize(sel);

    std::vector<int> num(system.configuration().num_particles());
    for (int i = 0; i < 200; ++i) {
      sel2.sel(&system, ran.get());
      ++num[sel2.mobile().particle_index(0)];
    }
    // INFO(feasst_str(num));
    if (max_dist == 1.5) {
      // 100 100 0
      EXPECT_EQ(0, num[2]);
      EXPECT_GT(num[0], 50);
      EXPECT_GT(num[1], 50);
    } else {
      // 133, 33, 34
      EXPECT_GT(num[0], 100);
      EXPECT_LT(2*num[1], num[0]);
      EXPECT_LT(2*num[2], num[0]);
    }
  }
}

}  // namespace feasst
