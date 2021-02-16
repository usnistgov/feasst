#include "utils/test/utils.h"
#include "cluster/include/trial_avb4.h"

namespace feasst {

TEST(TrialAVB4, serialize) {
  auto trial = MakeTrialAVB4({{"neighbor_index", "0"}});
  Trial trial2 = test_serialize(*trial);
}

}  // namespace feasst
