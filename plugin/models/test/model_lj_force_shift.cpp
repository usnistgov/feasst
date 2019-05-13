#include "utils/test/utils.h"
#include "models/include/model_lj_force_shift.h"
#include "configuration/include/configuration.h"

namespace feasst {

TEST(ModelLJForceShift, analytical) {
  Configuration config;
  config.add_particle_type("../forcefield/data.lj");
  EXPECT_NEAR(3., config.model_params().mixed_cutoff()[0][0], NEAR_ZERO);
  auto force_shift = std::make_shared<ModelLJForceShift>();
  force_shift->precompute(config.model_params());
  EXPECT_NEAR(-0, force_shift->energy(3*3, 0, 0, config.model_params()), NEAR_ZERO);
  EXPECT_NEAR(-0.005365534353215080, force_shift->energy(2.5*2.5, 0, 0, config.model_params()), NEAR_ZERO);
}

TEST(ModelLJForceShift, serialize) {
  Configuration config;
  config.add_particle_type("../forcefield/data.spce");
  auto shift = std::make_shared<ModelLJForceShift>();
  shift->set_alpha(12);
  shift->set_hard_sphere_threshold(0.3);
  shift->precompute(config.model_params());
  std::shared_ptr<Model> model2 = test_serialize<ModelLJForceShift, Model>(*shift,
    "ModelLJForceShift 763 0.089999999999999997 713 12 923 1 ModelParam 1 2 0 0 2 2 -2.6332331818264547e-06 -0 2 -0 -0 1 1 ModelParam 1 2 0 0 2 2 3.1598766187506022e-06 0 2 0 0 1 1 ");
}

}  // namespace feasst
