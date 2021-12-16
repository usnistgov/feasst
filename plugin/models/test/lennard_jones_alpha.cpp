#include "utils/test/utils.h"
#include "models/include/lennard_jones_force_shift.h"
#include "configuration/include/configuration.h"

namespace feasst {

TEST(LennardJonesAlpha, analytical) {
  Configuration config;
  config.add_particle_type("../forcefield/lj.fstprt");
  auto model1 = std::make_shared<LennardJones>();
  auto model2 = std::make_shared<LennardJonesAlpha>();
  model1->precompute(config.model_params());
  model2->precompute(config.model_params());
  EXPECT_NEAR(model1->energy(3.*3., 0, 0, config.model_params()),
              model2->energy(3.*3., 0, 0, config.model_params()), NEAR_ZERO);
}

double en(const double distance) {
  const double tmp = std::pow(1./(distance+0.5), 6);
  return 4*tmp*(tmp - 1.);
}

TEST(LennardJonesAlpha, analytical_ref) {
  Configuration config;
  config.add(MakeDeltaSigma());
  config.add_particle_type("../plugin/models/forcefield/ljdelta.fstprt");
  auto model = MakeLennardJonesAlpha();
  model->precompute(config.model_params());
  EXPECT_NEAR(en(3), model->energy(3.*3., 0, 0, config.model_params()), NEAR_ZERO);
  EXPECT_NEAR(en(1.2), -0.15885125916246515, NEAR_ZERO);
  EXPECT_NEAR(en(1.2), model->energy(1.2*1.2, 0, 0, config.model_params()), NEAR_ZERO);
}

TEST(LennardJonesAlpha, serialize) {
  Configuration config;
  config.add_particle_type("../forcefield/lj.fstprt");
  auto model = MakeLennardJonesAlpha({{"alpha", "12"},
                                      {"hard_sphere_threshold", "0.3"}});
  model->precompute(config.model_params());
  std::shared_ptr<Model> model2 = test_serialize<LennardJonesAlpha, Model>(*model, "LennardJonesAlpha 2094 0 1 2 3 763 0.089999999999999997 713 12 -1 0 -1 1.0594630943592953 ");
}

TEST(LennardJonesAlpha, analytical_lambda) {
  Configuration config;
  config.add(MakeLambda());
  config.add(MakeDeltaSigma());
  config.add_particle_type("../plugin/models/forcefield/ljlambda.fstprt");
  auto model = MakeLennardJonesAlpha({{"lambda", "true"}});
  model->precompute(config.model_params());
  //EXPECT_NEAR(en(3), model->energy(3.*3., 0, 0, config.model_params()), NEAR_ZERO);
  EXPECT_NEAR(-0.002739720872119390, model->energy(2.5*2.5, 0, 0, config.model_params()), NEAR_ZERO);
  EXPECT_NEAR(-0.001087390195827500, model->energy(3*3, 0, 0, config.model_params()), NEAR_ZERO);
}

}  // namespace feasst
