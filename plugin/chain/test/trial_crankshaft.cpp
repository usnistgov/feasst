#include "utils/test/utils.h"
#include "chain/include/trial_crankshaft.h"

namespace feasst {

TEST(TrialCrankshaft, serialize) {
  TrialCrankshaft add;
  TrialCrankshaft add2 = test_serialize(add);
}

}  // namespace feasst
