#include "utils/test/utils.h"
#include "chain/include/perturb_reptate.h"

namespace feasst {

TEST(PerturbReptate, serialize) {
  PerturbReptate add;
  PerturbReptate add2 = test_serialize(add);
}

}  // namespace feasst
