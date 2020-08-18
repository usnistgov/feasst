#include "utils/test/utils.h"
#include "configuration/include/configuration.h"
#include "ewald/include/utils.h"

namespace feasst {

TEST(Configuration, synchronize) {
  System s1 = spce();
  s1.get_configuration()->add_particle_of_type(0);
  s1.energy();
  EXPECT_NEAR(s1.configuration().particle(0).site(0).property("eikrx2"), 1., NEAR_ZERO);
  System s2 = test_serialize(s1);
  Select part(0, s1.configuration().particle(0));
  Position disp({0.5, 0.5, 0.5});
  s1.get_configuration()->displace_particle(part, disp);
  s1.energy();
  EXPECT_NEAR(s1.configuration().particle(0).site(0).position().coord(0), 0.5, NEAR_ZERO);
  EXPECT_NEAR(s1.configuration().particle(0).site(0).property("eikrx2"), 0.95105651629515364, NEAR_ZERO);
  EXPECT_NEAR(s2.configuration().particle(0).site(0).position().coord(0), 0., NEAR_ZERO);
  EXPECT_NEAR(s2.configuration().particle(0).site(0).property("eikrx2"), 1, NEAR_ZERO);
  s2.synchronize_(s1, part);
  EXPECT_NEAR(s2.configuration().particle(0).site(0).position().coord(0), 0.5, NEAR_ZERO);
  EXPECT_NEAR(s2.configuration().particle(0).site(0).property("eikrx2"), 0.95105651629515364, NEAR_ZERO);
}

}  // namespace feasst
