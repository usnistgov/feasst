#include <gtest/gtest.h>
#include "core/include/model_lj.h"
#include "core/include/configuration.h"

namespace feasst {

TEST(ModelLJ, analytical) {
  Configuration config;
  config.add_particle_type("../forcefield/data.lj");
  auto model = std::make_shared<ModelLJ>();
  EXPECT_NEAR(-0.005479441744238780, model->energy(3.*3., 0, 0, config.model_params()), NEAR_ZERO);

  // test lj alpha
  auto model1 = std::make_shared<ModelLJ>();
  auto model2 = std::make_shared<ModelLJAlpha>();
  EXPECT_NEAR(model->energy(3.*3., 0, 0, config.model_params()),
              model->energy(3.*3., 0, 0, config.model_params()), NEAR_ZERO);

  // test lj cut shift
  auto shift = std::make_shared<ModelLJCutShift>();
  shift->precompute(config.model_params());
  EXPECT_NEAR(0., shift->energy(3*3, 0, 0, config.model_params()), NEAR_ZERO);
  EXPECT_NEAR(-0.010837449391761200, shift->energy(2.5*2.5, 0, 0, config.model_params()), NEAR_ZERO);

  // test lj force shift
  EXPECT_NEAR(3., config.model_params().mixed_cutoff()[0][0], NEAR_ZERO);
  auto force_shift = std::make_shared<ModelLJForceShift>();
  force_shift->precompute(config.model_params());
  EXPECT_NEAR(-0, force_shift->energy(3*3, 0, 0, config.model_params()), NEAR_ZERO);
  EXPECT_NEAR(-0.005365534353215080, force_shift->energy(2.5*2.5, 0, 0, config.model_params()), NEAR_ZERO);

  // test wca
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


}  // namespace feasst
