#include "utils/test/utils.h"
#include "chain/include/ghost_trial_grow.h"

namespace feasst {

TEST(GhostTrialGrow, serialize) {
  GhostTrialGrow reject;
  GhostTrialGrow reject2 = test_serialize(reject);
}

}  // namespace feasst
