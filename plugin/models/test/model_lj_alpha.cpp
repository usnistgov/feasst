#include <gtest/gtest.h>
#include "models/include/model_lj_force_shift.h"
#include "core/include/configuration.h"

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
  std::stringstream ss, ss2;
  model.serialize(ss);
  EXPECT_EQ("ModelLJAlpha 763 0.089999999999999997 713 12 ", ss.str());
  std::shared_ptr<Model> model2 = ModelLJAlpha().deserialize(ss);
  model2->serialize(ss2);
  EXPECT_EQ(ss.str(), ss2.str());
}

}  // namespace feasst
