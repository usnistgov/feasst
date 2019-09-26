#include "utils/test/utils.h"
#include "monte_carlo/include/perturb_distance.h"

namespace feasst {

TEST(PerturbDistance, serialize) {
  PerturbDistance add;
  PerturbDistance add2 = test_serialize(add);
}

}  // namespace feasst
