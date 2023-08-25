#include "utils/test/utils.h"
#include "chain/include/perturb_position_swap.h"

namespace feasst {

TEST(PerturbPositionSwap, serialize) {
  PerturbPositionSwap perturb;
  PerturbPositionSwap perturb2 = test_serialize(perturb);
}

}  // namespace feasst
