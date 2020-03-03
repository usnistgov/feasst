#include "utils/test/utils.h"
#include "cluster/include/trial_select_cluster.h"
#include "system/test/system_test.h"
#include "math/include/random_mt19937.h"
#include "monte_carlo/include/perturb_translate.h"
#include "cluster/include/energy_map_all.h"
#include "cluster/include/perturb_rotate_com.h"

namespace feasst {

inline double en_lj(const double dist) {
  return 4*(pow(dist, -12) - pow(dist, -6));
}

TEST(TrialSelectCluster, serialize) {
  //MonteCarlo mc;
  //mc_lj(&mc);
  System sys;
  {
    Configuration config(MakeDomain({{"cubic_box_length", "8"}}),
                         {{"particle_type", "../forcefield/data.lj"}});
    config.add_particle_of_type(0);
    config.add_particle_of_type(0);
    config.add_particle_of_type(0);
    config.update_positions({{0, 0, 0},
                             {-1.25, 0, 0},
                             {3, 3, 3}});
    sys.add(config);
  }
  sys.add(Potential(MakeLennardJones(),
                    MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
  sys.energy();
  sys.finalize();
  TrialSelectCluster add(MakeClusterCriteria());

  // select all clusters
  { std::vector<SelectList> clusters = add.select_clusters(sys);
    EXPECT_EQ(2, static_cast<int>(clusters.size()));
    EXPECT_EQ(2, clusters[0].num_particles());
    EXPECT_EQ(1, clusters[1].num_particles());
    EXPECT_EQ(2, clusters[1].particle_index(0));
  }

  auto ran = MakeRandomMT19937();
  ran = MakeRandomMT19937({{"seed", "1580154124"}});
  add.select(Select(), &sys, ran.get());
  // DEBUG(add.mobile().str());

  if (add.mobile().num_particles() == 1) {
    EXPECT_EQ(2, add.mobile().particle_index(0));
  } else {
    EXPECT_EQ(2, add.mobile().num_particles());
    if (add.mobile().particle_index(0) == 0) {
      EXPECT_EQ(1, add.mobile().particle_index(1));
    } else {
      EXPECT_EQ(1, add.mobile().particle_index(0));
    }
    DEBUG(add.mobile().str());
    EXPECT_NEAR(en_lj(1.25), sys.perturbed_energy(add.mobile()), NEAR_ZERO);
  }

  TrialSelectCluster add2 = test_serialize(add);

  PerturbTranslate trans;
  trans.perturb(&sys, &add, ran.get());
  if (add.mobile().num_particles() == 1) {
    DEBUG(sys.configuration().particle(2).site(0).position().str());
    EXPECT_NEAR(0, sys.perturbed_energy(add.mobile()), NEAR_ZERO);
  } else {
    DEBUG(sys.configuration().particle(0).site(0).position().str());
    DEBUG(sys.configuration().particle(1).site(0).position().str());
    EXPECT_NEAR(en_lj(1.25), sys.perturbed_energy(add.mobile()), NEAR_ZERO);
  }

  auto rotate = MakePerturbRotateCOM({{"tunable_param", "50"}});
  rotate->perturb(&sys, &add, ran.get());
  if (add.mobile().num_particles() == 1) {
    DEBUG(sys.configuration().particle(2).site(0).position().str());
    EXPECT_NEAR(0, sys.perturbed_energy(add.mobile()), NEAR_ZERO);
  } else {
    DEBUG(sys.configuration().particle(0).site(0).position().str());
    DEBUG(sys.configuration().particle(1).site(0).position().str());
    EXPECT_NEAR(en_lj(1.25), sys.perturbed_energy(add.mobile()), NEAR_ZERO);
  }
}

}  // namespace feasst
