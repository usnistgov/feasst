#include <gtest/gtest.h>
#include "models/include/model_lj_cut_shift.h"
#include "core/include/configuration.h"

namespace feasst {

TEST(ModelLJCutShift, analytical) {
  Configuration config;
  config.add_particle_type("../forcefield/data.lj");
  auto shift = std::make_shared<ModelLJCutShift>();
  shift->precompute(config.model_params());
  EXPECT_NEAR(0., shift->energy(3*3, 0, 0, config.model_params()), NEAR_ZERO);
  EXPECT_NEAR(-0.010837449391761200, shift->energy(2.5*2.5, 0, 0, config.model_params()), NEAR_ZERO);
}

TEST(ModelLJCutShift, wca) {
  Configuration config;
  config.add_particle_type("../forcefield/data.lj");
  auto wca = std::make_shared<ModelLJCutShift>();
  ModelParams wca_params = config.model_params();
  EXPECT_NEAR(3., config.model_params().mixed_cutoff()[0][0], NEAR_ZERO);
  wca->set_wca(0, 0, &wca_params);
  EXPECT_NEAR(3., config.model_params().mixed_cutoff()[0][0], NEAR_ZERO);
  wca->precompute(wca_params);
  const double r_wca = pow(2, 1./6.);
  EXPECT_NEAR(0., wca->energy(r_wca*r_wca, 0, 0, wca_params), NEAR_ZERO);
  EXPECT_NEAR(1., wca->energy(1., 0, 0, wca_params), NEAR_ZERO);
}

TEST(ModelLJCutShift, serialize) {
  Configuration config;
  config.add_particle_type("../forcefield/data.spce");
  auto shift = std::make_shared<ModelLJCutShift>();
  shift->set_alpha(12);
  shift->set_hard_sphere_threshold(0.3);
  shift->precompute(config.model_params());
  std::stringstream ss, ss2;
  shift->serialize(ss);
  std::string expected("ModelLJCutShift 763 0.089999999999999997 713 12 644 1 ModelParam 1 2 0 0 2 2 -2.6332331818264547e-06 -0 2 -0 -0 1 1 ");
  EXPECT_EQ(expected, ss.str());
  std::shared_ptr<Model> model2 = ModelLJAlpha().deserialize(ss);
  model2->serialize(ss2);
  EXPECT_EQ(expected, ss2.str());
}

}  // namespace feasst
