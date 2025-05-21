#include "utils/test/utils.h"
#include "monte_carlo/test/monte_carlo_utils.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/configuration.h"
#include "configuration/include/domain.h"
#include "system/include/system.h"
#include "system/include/visit_model.h"
#include "system/include/visit_model_inner.h"
#include "system/include/potential.h"
#include "system/include/lennard_jones.h"
#include "cluster/include/energy_map_all.h"
#include "cluster/include/select_particle_avb.h"
#include "cluster/include/perturb_add_avb.h"

namespace feasst {

TEST(PerturbAddAVB, gce_add) {
  System system;
  {
    auto config = MakeConfiguration({{"cubic_side_length", "8"},
      {"particle_type", "../particle/lj.txt"},
      {"add_particles_of_type0", "1"}});
    config->update_positions({{0, 0, 0}});
    system.add(config);
  }
  system.add(MakePotential(MakeLennardJones(),
                       MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
  system.energy();
  system.finalize();

  auto ran = MakeRandomMT19937();
  //ran = MakeRandomMT19937({{"seed", "1580154124"}});

  const double max_dist = 0.1;
  system.add(MakeNeighborCriteria({{"maximum_distance", str(max_dist)}}));

  SelectParticleAVB sel({
    {"neighbor_index", "0"},
    {"grand_canonical", "true"},
    {"particle_type", "0"}});
  sel.precompute(&system);
  auto sel2 = test_serialize(sel);

  auto add = MakePerturbAddAVB({{"neighbor_index", "0"}});
  add->precompute(&sel2, &system);
  auto add2 = test_serialize(*add);

  EXPECT_TRUE(sel2.is_ghost());
  sel2.sel(&system, ran.get());
  add2.perturb(&system, &sel2, ran.get());
  add2.finalize(&system);

  EXPECT_EQ(sel2.mobile().particle_index(0), 1);
  EXPECT_EQ(sel2.mobile().site_index(0, 0), 0);
  EXPECT_EQ(sel2.anchor().particle_index(0), 0);
  EXPECT_EQ(sel2.anchor().site_index(0, 0), 0);
  EXPECT_EQ(system.configuration().num_particles(), 2);
  EXPECT_LT(system.configuration().particle(1).site(0).position().distance(),
            max_dist);
}

}  // namespace feasst
