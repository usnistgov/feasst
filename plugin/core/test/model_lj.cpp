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

}  // namespace feasst
