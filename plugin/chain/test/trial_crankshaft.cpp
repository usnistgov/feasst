#include "utils/test/utils.h"
#include "chain/include/trials.h"

namespace feasst {

TEST(TrialCrankshaft, serialize) {
  auto trial = MakeTrialCrankshaft();
  Trial trial2 = test_serialize(*trial);
}

}  // namespace feasst
