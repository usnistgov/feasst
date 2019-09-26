#include "utils/test/utils.h"
#include "chain/include/perturb_pivot.h"

namespace feasst {

TEST(PerturbPivot, serialize) {
  PerturbPivot add;
  PerturbPivot add2 = test_serialize(add);
}

}  // namespace feasst
