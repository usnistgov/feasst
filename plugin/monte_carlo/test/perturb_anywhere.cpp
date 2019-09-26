#include "utils/test/utils.h"
#include "monte_carlo/include/perturb_anywhere.h"

namespace feasst {

TEST(PerturbAnywhere, serialize) {
  PerturbAnywhere add;
  PerturbAnywhere add2 = test_serialize(add);
}

}  // namespace feasst
