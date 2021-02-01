#include <sstream>
#include "utils/test/utils.h"
#include "system/include/lennard_jones.h"
#include "configuration/include/configuration.h"

namespace feasst {

TEST(LennardJones, analytical) {
  Configuration config;
  config.add_particle_type("../forcefield/data.lj");
  auto model = std::make_shared<LennardJones>();
  EXPECT_NEAR(-0.005479441744238780, model->energy(3.*3., 0, 0, config.model_params()), NEAR_ZERO);
}

TEST(LennardJones, serialize) {
  auto model = MakeLennardJones({{"hard_sphere_threshold", "0.3"}});
  auto model2 = test_serialize<LennardJones, Model>(*model,
    "LennardJones 763 0.089999999999999997 ");
}

}  // namespace feasst
