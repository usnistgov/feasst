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
}

}  // namespace feasst
