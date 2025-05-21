#include <sstream>
#include "utils/test/utils.h"
#include "system/include/lennard_jones.h"
#include "configuration/include/configuration.h"

namespace feasst {

TEST(LennardJones, analytical) {
  Configuration config;
  config.add_particle_type("../particle/lj.txt");
  auto model = std::make_shared<LennardJones>();
  model->precompute(config.model_params());
  EXPECT_NEAR(-0.005479441744238780, model->energy(3.*3., 0, 0, config.model_params()), NEAR_ZERO);
}

TEST(LennardJones, serialize) {
  Configuration config;
  config.add_particle_type("../particle/lj.txt");
  auto model = MakeLennardJones({{"hard_sphere_threshold", "0.3"}});
  model->precompute(config.model_params());
  auto model2 = test_serialize<LennardJones, Model>(*model,
    "LennardJones 2094 1 0 2 -1 763 0.089999999999999997 ");
}

}  // namespace feasst
