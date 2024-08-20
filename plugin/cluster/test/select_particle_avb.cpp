#include "utils/test/utils.h"
#include "monte_carlo/test/monte_carlo_utils.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"
#include "system/include/visit_model_cell.h"
#include "system/include/visit_model_inner.h"
#include "system/include/potential.h"
#include "system/include/thermo_params.h"
#include "system/include/lennard_jones.h"
#include "monte_carlo/include/perturb_translate.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/energy_map_all.h"
#include "cluster/include/perturb_rotate_com.h"

namespace feasst {

TEST(SelectParticleAVB, serialize) {
  //for (const double max_dist : {1.5}) {
  for (const double max_dist : {1.5, 3.}) {
  //for (const double max_dist : {3.}) {
    System system;
    {
      auto config = MakeConfiguration({{"cubic_side_length", "8"},
        {"particle_type", "../particle/lj.fstprt"},
        {"add_particles_of_type0", "3"}});
      config->update_positions({{0, 0, 0},
                                {-1.25, 0, 0},
                                {2.9, 0, 0}});
      system.add(config);
    }
    system.add(MakePotential(MakeLennardJones(),
      MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
    system.add(MakeNeighborCriteria({{"maximum_distance", str(max_dist)}}));
    system.energy();
    system.finalize();

    auto ran = MakeRandomMT19937();
  //ran = MakeRandomMT19937({{"seed", "1580154124"}});

    SelectParticleAVB sel(
      {{"grand_canonical", "true"},
       {"neighbor_index", "0"},
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
