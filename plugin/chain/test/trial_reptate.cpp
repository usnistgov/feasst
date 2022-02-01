#include "utils/test/utils.h"
#include "chain/include/trial_reptate.h"

namespace feasst {

TEST(TrialReptate, serialize) {
  auto trial = MakeTrialReptate({{"max_length", "1"}});
  Trial trial2 = test_serialize(*trial);
}

}  // namespace feasst
