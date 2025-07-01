#include "utils/test/utils.h"
#include "configuration/include/model_params.h"
#include "configuration/include/configuration.h"

namespace feasst {

TEST(Configuration, parameter_not_in_fstprt) {
  auto config = MakeConfiguration({{"particle_type", "lj:../particle/lj_new.txt"},
                                   {"patch_angle", "3"}});
  EXPECT_NEAR(config->model_params().select("patch_angle").value(0), 3, NEAR_ZERO);
}


}  // namespace feasst
