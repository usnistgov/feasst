#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "system/include/lennard_jones.h"
#include "cluster/include/energy_map_all.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/perturb_add_avb.h"

namespace feasst {

TEST(PerturbAddAVB, gce_add) {
  System system;
  {
    Configuration config(MakeDomain({{"cubic_box_length", "8"}}),
                         {{"particle_type", "../forcefield/data.lj"}});
    config.add_particle_of_type(0);
    config.update_positions({{0, 0, 0}});
    system.add(config);
  }
  system.add(Potential(MakeLennardJones(),
                       MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
  system.energy();
  system.finalize();

  auto ran = MakeRandomMT19937();
  //ran = MakeRandomMT19937({{"seed", "1580154124"}});

  const double max_dist = 0.1;
  auto neigh_crit = MakeNeighborCriteria({{"maximum_distance", str(max_dist)}});

  SelectParticleAVB sel(
    neigh_crit,
    {{"grand_canonical", "true"},
     {"particle_type", "0"}});
  sel.precompute(&system);
  auto sel2 = test_serialize(sel);

  PerturbAddAVB add(neigh_crit);
  add.precompute(&sel2, &system);
  auto add2 = test_serialize(add);

  EXPECT_TRUE(sel2.is_ghost());
  sel2.sel(&system, ran.get());
  add2.perturb(&system, &sel2, ran.get());

  EXPECT_EQ(sel2.mobile().particle_index(0), 1);
  EXPECT_EQ(sel2.mobile().site_index(0, 0), 0);
  EXPECT_EQ(sel2.anchor().particle_index(0), 0);
  EXPECT_EQ(sel2.anchor().site_index(0, 0), 0);
  EXPECT_EQ(system.configuration().num_particles(), 2);
  EXPECT_LT(system.configuration().particle(1).site(0).position().distance(),
            max_dist);
}

}  // namespace feasst
