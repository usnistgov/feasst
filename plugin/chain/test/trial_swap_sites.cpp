#include "utils/test/utils.h"
#include "configuration/include/domain.h"
#include "system/include/lennard_jones.h"
#include "math/include/random_mt19937.h"
#include "monte_carlo/include/metropolis.h"
#include "chain/include/trials.h"

namespace feasst {

TEST(TrialSwapSites, swap) {
  System sys;
  {
    Configuration config(MakeDomain({{"cubic_box_length", "20"}}),
      {{"particle_type", "../forcefield/chain10_3types.fstprt"}});
    config.add_particle_of_type(0);
    sys.add(config);
  }
  sys.add(MakePotential(MakeLennardJones()));
  sys.energy();
  sys.set(MakeThermoParams({{"beta", "1.0"}}));
  auto random = MakeRandomMT19937();
  auto criteria = MakeMetropolis();
  auto trial = MakeTrialSwapSites({{"particle_type", "0"}, {"site_type1", "0"}, {"site_type2", "1"}});
  trial->attempt(criteria.get(), &sys, random.get());
  EXPECT_EQ(trial->num_success(), 1);
  EXPECT_EQ(sys.configuration().particle(0).site(0).type(), 0);
  const int site_index = trial->stage(0).trial_select().mobile().site_index(0, 0);
  EXPECT_TRUE(site_index == 1 ||
              site_index == 3 ||
              site_index == 4 ||
              site_index == 5 ||
              site_index == 6 ||
              site_index == 7 ||
              site_index == 8 ||
              site_index == 9);
  EXPECT_EQ(sys.configuration().particle(0).site(site_index).type(), 1);
}

}  // namespace feasst
