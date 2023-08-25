#include "utils/test/utils.h"
#include "chain/include/perturb_to_anchor.h"

namespace feasst {

TEST(PerturbToAnchor, serialize) {
  PerturbToAnchor perturb;
  PerturbToAnchor perturb2 = test_serialize(perturb);
}

}  // namespace feasst
