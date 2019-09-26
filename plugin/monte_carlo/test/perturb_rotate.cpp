#include "utils/test/utils.h"
#include "monte_carlo/include/perturb_rotate.h"

namespace feasst {

TEST(PerturbRotate, serialize) {
  PerturbRotate add;
  PerturbRotate add2 = test_serialize(add);
}

}  // namespace feasst
