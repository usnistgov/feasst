#include "utils/test/utils.h"
#include "morph/include/trial_position_swap.h"

namespace feasst {

TEST(TrialPositionSwap, serialize) {
  auto obj = std::make_shared<TrialPositionSwap>(argtype({{"particle_type", "0"}, {"particle_type_morph", "1"}}));
  std::shared_ptr<Trial> trial2 = test_serialize<TrialPositionSwap, Trial>(*obj);
}

}  // namespace feasst
