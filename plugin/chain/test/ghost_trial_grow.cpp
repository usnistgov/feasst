#include "utils/test/utils.h"
#include "chain/include/ghost_trial_grow.h"

namespace feasst {

TEST(GhostTrialGrow, serialize) {
  GhostTrialGrow reject;
  auto reject2 = test_serialize_unique(reject);
}

}  // namespace feasst
