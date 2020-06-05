#include "utils/test/utils.h"
#include "system/include/utils.h"

namespace feasst {

TEST(System, utils) {
  System lj = lennard_jones({{"dual_cut", "1."}});
  EXPECT_FALSE(lj.potential(0).are_model_params_overridden());
  EXPECT_NEAR(3,
    lj.potential(0).model_params(lj.configuration()).cutoff().value(0),
    NEAR_ZERO);
  EXPECT_TRUE(lj.reference(0, 0).are_model_params_overridden());
  EXPECT_NEAR(1., lj.reference(0, 0).model_params().cutoff().value(0),
    NEAR_ZERO);
}

}  // namespace feasst
