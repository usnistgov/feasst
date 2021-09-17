#include "utils/test/utils.h"
#include "models/include/lennard_jones_force_shift.h"
#include "configuration/include/configuration.h"

namespace feasst {

TEST(LennardJonesAlpha, analytical) {
  Configuration config;
  config.add_particle_type("../forcefield/lj.fstprt");
  auto model1 = std::make_shared<LennardJones>();
  auto model2 = std::make_shared<LennardJonesAlpha>();
  EXPECT_NEAR(model1->energy(3.*3., 0, 0, config.model_params()),
              model2->energy(3.*3., 0, 0, config.model_params()), NEAR_ZERO);
}

TEST(LennardJonesAlpha, serialize) {
  auto model = MakeLennardJonesAlpha({{"alpha", "12"},
                                      {"hard_sphere_threshold", "0.3"}});
  std::shared_ptr<Model> model2 = test_serialize<LennardJonesAlpha, Model>(*model, "LennardJonesAlpha 763 0.089999999999999997 713 12 ");
}

}  // namespace feasst
