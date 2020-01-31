#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "system/include/energy_map_all.h"
#include "system/include/lennard_jones.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/metropolis.h"
#include "steppers/include/log.h"
#include "steppers/include/movie.h"
#include "steppers/include/check_energy.h"
#include "cluster/include/utils_cluster.h"
#include "cluster/include/trial_select_cluster.h"

namespace feasst {

/// Without single particle translations, rigid cluster moves should reject
/// cluster coalescence and breakup to satisfy detailed balance.
TEST(MonteCarlo, cluster) {
  for (bool single_particle_translate : {true, false}) {
    MonteCarlo monte_carlo;
    // monte_carlo.set(MakeRandomMT19937({{"seed", "default"}}));
    monte_carlo.add(Configuration({{"cubic_box_length", "8"},
                                   {"particle_type", "../forcefield/data.lj"}}));
    monte_carlo.add(Potential(MakeLennardJones(),
                              MakeVisitModel(MakeVisitModelInner(MakeEnergyMapAll()))));
    monte_carlo.add(MakeMetropolis({{"beta", "40"}, {"chemical_potential", "1."}}));
    if (single_particle_translate) monte_carlo.add(MakeTrialTranslate());
    monte_carlo.seek_num_particles(3);
    auto cluster_criteria = MakeClusterCriteria({{"energy_maximum", "-0.5"}});
    TrialSelectCluster scluster(cluster_criteria);
    scluster.select_cluster(0, monte_carlo.system());
    const int cluster_size = scluster.mobile().num_particles();
    add_rigid_cluster_trials(&monte_carlo,
      cluster_criteria,
      {{"tunable_param", "50"}});
    const int steps_per = 1e3;
    monte_carlo.add(MakeMovie({{"steps_per", str(steps_per)},
                               {"file_name", "cluster.xyz"},
                               {"clear_file", "true"}}));
    monte_carlo.add(MakeLog({{"steps_per", str(steps_per)},
                             {"file_name", "cluster.txt"},
                             {"clear_file", "true"}}));
    monte_carlo.add(MakeCheckEnergy({{"steps_per", str(steps_per)},
                                     {"tolerance", str(1e-8)}}));
    monte_carlo.attempt(1e3);
    scluster.select_cluster(0, monte_carlo.system());
    if (single_particle_translate) {
      EXPECT_EQ(scluster.mobile().num_particles(),
                monte_carlo.configuration().num_particles());
    } else {
      EXPECT_EQ(scluster.mobile().num_particles(), cluster_size);
    }
  }
}

  //  // test rotation through PBCs
  //  monte_carlo.get_system()->get_configuration()->update_positions({
  //    {4.0, 0.769800358919501000, 0},
  //    { 4.666666666666666666, -0.384900179459751000, 0},
  //    { 3.333333333333333334, -0.384900179459751000, 0}});
//    monte_carlo.get_system()->energy();
//    monte_carlo.add(MakeMetropolis({{"beta", "40"}, {"chemical_potential", "1."}}));
}  // namespace feasst
