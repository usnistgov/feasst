#include "utils/test/utils.h"
#include "rxvb/include/trial_swap_position.h"

namespace feasst {

TEST(TrialSwapPosition, serialize) {
  auto trial = std::make_unique<TrialSwapPosition>(argtype({{"particle_types", "1,2"}}));
  auto trial2 = test_serialize(trial);
}

}  // namespace feasst
