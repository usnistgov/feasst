#include "utils/test/utils.h"
#include "chain/include/perturb_crankshaft.h"

namespace feasst {

TEST(PerturbCrankshaft, serialize) {
  PerturbCrankshaft add;
  PerturbCrankshaft add2 = test_serialize(add);
}

}  // namespace feasst
