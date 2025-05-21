#include "utils/test/utils.h"
#include "monte_carlo/test/monte_carlo_utils.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/configuration.h"
#include "configuration/include/domain.h"
#include "system/include/system.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "chain/include/perturb_site_type.h"

namespace feasst {

TEST(PerturbSiteType, serialize) {
  System sys;
  sys.add(MakeConfiguration({{"cubic_side_length", "20"},
    {"particle_type", "../particle/chain10_3types.txt"},
    {"add_particles_of_type0", "1"}}));
  const Configuration& config = sys.configuration();
  auto morph = MakePerturbSiteType({{"type", "0"}});
  auto sel = MakeTrialSelectParticle({{"particle_type", "0"}, {"site", "1"}});
  sel->precompute(&sys);
  auto random = MakeRandomMT19937();
  sel->select(Select(), &sys, random.get(), NULL);
  EXPECT_EQ(config.particle(0).site(1).type(), 1);
  morph->set_site_type(&sys, *sel, 0);
  EXPECT_EQ(config.particle(0).site(1).type(), 0);
  PerturbSiteType morph2 = test_serialize(*morph);
}

}  // namespace feasst
