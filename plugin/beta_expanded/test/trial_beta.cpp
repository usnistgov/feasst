#include "utils/test/utils.h"
#include "beta_expanded/include/trial_beta.h"

namespace feasst {

TEST(TrialBeta, serialize) {
  auto trial = MakeTrialBeta({{"fixed_beta_change", "1"}});
  Trial trial2 = test_serialize(*trial);
}

}  // namespace feasst
