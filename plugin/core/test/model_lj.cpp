#include <sstream>
#include <gtest/gtest.h>
#include "core/include/model_lj.h"
#include "core/include/configuration.h"

namespace feasst {

TEST(ModelLJ, analytical) {
  Configuration config;
  config.add_particle_type("../forcefield/data.lj");
  auto model = std::make_shared<ModelLJ>();
  EXPECT_NEAR(-0.005479441744238780, model->energy(3.*3., 0, 0, config.model_params()), NEAR_ZERO);
}

TEST(ModelLJ, serialize) {
  ModelLJ model;
  model.set_hard_sphere_threshold(0.3);
  std::stringstream ss, ss2;
  model.serialize(ss);
  EXPECT_EQ("ModelLJ 763 0.089999999999999997 ", ss.str());
  std::shared_ptr<Model> model2 = ModelLJ().deserialize(ss);
  model2->serialize(ss2);
  EXPECT_EQ(ss.str(), ss2.str());
}

}  // namespace feasst
