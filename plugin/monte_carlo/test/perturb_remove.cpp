#include "utils/test/utils.h"
#include "monte_carlo/include/perturb_remove.h"

namespace feasst {

TEST(PerturbRemove, serialize) {
  PerturbRemove add;
  PerturbRemove add2 = test_serialize(add);
}

}  // namespace feasst
