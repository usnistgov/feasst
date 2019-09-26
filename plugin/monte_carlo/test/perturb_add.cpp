#include "utils/test/utils.h"
#include "monte_carlo/include/perturb_add.h"

namespace feasst {

TEST(PerturbAdd, serialize) {
  PerturbAdd add;
  PerturbAdd add2 = test_serialize(add);
}

}  // namespace feasst
