#include "utils/test/utils.h"
#include "chain/include/trial_reptate.h"

namespace feasst {

TEST(TrialReptate, serialize) {
  TrialReptate add;
  TrialReptate add2 = test_serialize(add);
}

}  // namespace feasst
