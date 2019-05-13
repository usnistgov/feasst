#include <sstream>
#include "utils/test/utils.h"
#include "system/include/model_lj.h"
#include "configuration/include/configuration.h"

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
  auto model2 = test_serialize<ModelLJ, Model>(model,
    "ModelLJ 763 0.089999999999999997 ");
}

}  // namespace feasst
