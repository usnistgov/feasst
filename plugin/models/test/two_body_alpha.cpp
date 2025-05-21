#include "utils/test/utils.h"
#include "configuration/include/configuration.h"
#include "models/include/two_body_alpha.h"

namespace feasst {

TEST(TwoBodyAlpha, serialize) {
  auto config = MakeConfiguration({
    {"particle_type0", "../particle/lj.txt"},
    {"sigma0", "2.5"},
    {"epsilon0", "3.7"}});
  auto model = MakeTwoBodyAlpha({
    {"alpha0", "12"}, {"s0", "1"},
    {"alpha1", "6"}, {"s1", "-1"}});
  model->precompute(config->model_params());
  EXPECT_NEAR(39.727706462144852,
    model->energy(2.*2., 0, 0, config->model_params()),
    NEAR_ZERO);
}

}  // namespace feasst
