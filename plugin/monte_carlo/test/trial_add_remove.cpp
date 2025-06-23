#include "utils/test/utils.h"
#include "monte_carlo/include/trial_add_remove.h"

namespace feasst {

TEST(TrialAddRemove, serialize) {
  auto trial = std::make_unique<TrialAddRemove>();
  auto trial2 = test_serialize_unique(*trial);
}

}  // namespace feasst
