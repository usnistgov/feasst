#include "utils/test/utils.h"
#include "monte_carlo/include/trials.h"

namespace feasst {

TEST(TrialAdd, serialize) {
  auto add = MakeTrialAdd();
  //Trial add2 = test_serialize(*add);
  std::shared_ptr<Trial> trial2 = test_serialize<TrialAdd, Trial>(*add);
  TRY(
    auto add3 = MakeTrialAdd({{"delay_add", "false"}});
    CATCH_PHRASE("TrialComputeAdd assumes delay_add is true");
  );
}

}  // namespace feasst
