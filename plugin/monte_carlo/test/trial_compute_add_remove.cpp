#include "utils/test/utils.h"
#include "monte_carlo/include/trial_compute_add_remove.h"

namespace feasst {

TEST(TrialComputeAddRemove, serialize) {
  auto compute = std::make_unique<TrialComputeAddRemove>();
  auto compute2 = test_serialize_unique(*compute);
}

}  // namespace feasst
