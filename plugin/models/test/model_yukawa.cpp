#include <gtest/gtest.h>
#include "models/include/model_yukawa.h"
#include "core/test/configuration_test.h"

namespace feasst {

TEST(ModelYukawa, analytical) {
  Configuration config = default_configuration();
  ModelYukawa model;
  model.set_kappa(2.);
  config.set_model_param("epsilon", 0, 0.5);
  EXPECT_NEAR(0.033833820809153176,
    model.energy(4., 0, 0, config.model_params()), NEAR_ZERO);

  // serialize
  std::stringstream ss;
  model.serialize(ss);
  std::shared_ptr<Model> model2 = ModelYukawa().deserialize(ss);
  std::stringstream ss2;
  model2->serialize(ss2);
  EXPECT_EQ(ss.str(), ss2.str());
}

}  // namespace feasst
