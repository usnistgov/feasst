#include "utils/test/utils.h"
#include "monte_carlo/include/trials.h"

namespace feasst {

TEST(TrialAdd, serialize) {
  auto add = MakeTrialAdd();
  Trial add2 = test_serialize(*add);
}

}  // namespace feasst
