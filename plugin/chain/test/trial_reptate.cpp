#include "utils/test/utils.h"
#include "chain/include/trial_reptate.h"

namespace feasst {

TEST(TrialReptate, serialize) {
  auto add = MakeTrialReptate({{"max_length", "1"}});
  TrialReptate add2 = test_serialize(*add);
}

}  // namespace feasst
