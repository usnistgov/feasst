#include "utils/test/utils.h"
#include "models/include/model_lj_force_shift.h"
#include "configuration/include/configuration.h"

namespace feasst {

TEST(ModelLJAlpha, analytical) {
  Configuration config;
  config.add_particle_type("../forcefield/data.lj");
  auto model1 = std::make_shared<ModelLJ>();
  auto model2 = std::make_shared<ModelLJAlpha>();
  EXPECT_NEAR(model1->energy(3.*3., 0, 0, config.model_params()),
              model2->energy(3.*3., 0, 0, config.model_params()), NEAR_ZERO);
}

TEST(ModelLJAlpha, serialize) {
  ModelLJAlpha model;
  model.set_alpha(12);
  model.set_hard_sphere_threshold(0.3);
  std::shared_ptr<Model> model2 = test_serialize<ModelLJAlpha, Model>(model, "ModelLJAlpha 763 0.089999999999999997 713 12 ");
}

}  // namespace feasst
