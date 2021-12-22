#include "utils/test/utils.h"
#include "models/include/yukawa.h"
#include "configuration/test/config_utils.h"

namespace feasst {

TEST(Yukawa, analytical) {
  Configuration config = two_particle_configuration();
  Yukawa model;
  model.precompute(config.model_params());
  model.set_kappa(2.);
  config.set_model_param("epsilon", 0, 0.5);
  EXPECT_NEAR(0.033833820809153176,
    model.energy(4., 0, 0, config.model_params()), NEAR_ZERO);
  std::shared_ptr<Model> model2 = test_serialize<Yukawa, Model>(model,
    "Yukawa 2094 1 0 2 -1 6505 2 ");
}

}  // namespace feasst
