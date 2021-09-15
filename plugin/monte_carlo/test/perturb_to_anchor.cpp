#include "utils/test/utils.h"
#include "monte_carlo/include/perturb_to_anchor.h"

namespace feasst {

TEST(PerturbToAnchor, serialize) {
  PerturbToAnchor perturb;
  PerturbToAnchor perturb2 = test_serialize(perturb);
}

}  // namespace feasst
