#include "utils/test/utils.h"
#include "rxvb/include/perturb_swap_position.h"

namespace feasst {

TEST(PerturbSwapPosition, serialize) {
  auto perturb = std::make_unique<PerturbSwapPosition>();
  auto perturb2 = test_serialize(perturb);
}

}  // namespace feasst
