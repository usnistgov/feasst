#include "utils/test/utils.h"
#include "monte_carlo/include/trial_stage.h"

namespace feasst {

TEST(TrialStage, serialize) {
  argtype args = {{"reference_index", "0"}};
  argtype tmp_args = args;
  auto stage = std::make_shared<TrialStage>(&tmp_args);
  TrialStage stage2 = test_serialize(*stage);
  argtype args2 = get_stage_args(&args);
  EXPECT_EQ("0", args2["reference_index"]);
}

}  // namespace feasst
