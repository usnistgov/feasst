#include "utils/test/utils.h"
#include "configuration/include/configuration.h"
#include "models/include/two_body_alpha.h"

namespace feasst {

TEST(TwoBodyAlpha, serialize) {
  auto config = MakeConfiguration({
    {"particle_type", "../particle/lj_new.txt"},
    {"sigmaLJ", "2.5"},
    {"epsilonLJ", "3.7"}});
  auto model = MakeTwoBodyAlpha({
    {"alpha", "12,6"}, {"s", "1,-1"}});
  model->precompute(*config);
  EXPECT_NEAR(39.727706462144852,
    model->energy(2.*2., 0, 0, config->model_params()),
    NEAR_ZERO);
}

}  // namespace feasst
