#include "utils/test/utils.h"
#include "configuration/include/domain.h"
#include "chain/include/perturb_site_type.h"

namespace feasst {

TEST(PerturbSiteType, serialize) {
  System sys;
  sys.add(*MakeConfiguration({{"cubic_side_length", "20"},
    {"particle_type", "../particle/chain10_3types.fstprt"},
    {"add_particles_of_type0", "1"}}));
  const Configuration& config = sys.configuration();
  auto morph = MakePerturbSiteType({{"type", "0"}});
  Select first_site;
  first_site.add_site(0, 1);
  EXPECT_EQ(config.particle(0).site(1).type(), 0);
  morph->set_site_type(&sys, first_site, 1);
  EXPECT_EQ(config.particle(0).site(1).type(), 1);
  PerturbSiteType morph2 = test_serialize(*morph);
}

}  // namespace feasst
