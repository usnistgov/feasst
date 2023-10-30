#include "utils/test/utils.h"
#include "cluster/include/trial_add_avb.h"

namespace feasst {

TEST(TrialAddAVB, serialize) {
  auto trial = MakeTrialAddAVB();
  Trial trial2 = test_serialize(*trial);
  TRY(
    auto trial3 = MakeTrialAddAVB({{"delay_add", "false"}});
    CATCH_PHRASE("ComputeAddAVB assumes delay_add is true");
  );
}

}  // namespace feasst
