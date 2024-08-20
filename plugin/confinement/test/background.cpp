#include "utils/test/utils.h"
#include "configuration/test/config_utils.h"
#include "system/include/system.h"
#include "system/include/lennard_jones.h"
#include "system/include/potential.h"
#include "confinement/include/background.h"

namespace feasst {

TEST(Background, constant) {
  Configuration config = lj_sample4();
  System system;
  system.add(std::make_shared<Configuration>(lj_sample4()));
  system.add(MakePotential(MakeLennardJones()));
  const double en_lj_expect = -16.790321304625856;
  const double constant = 10.;
  system.add(MakePotential(MakeBackground({{"constant", str(constant)}})));
  system.precompute();
  System system2 = test_serialize(system);
  const double en = system2.energy();
  INFO(en);
  EXPECT_DOUBLE_EQ(en, en_lj_expect + constant);
}

}  // namespace feasst
