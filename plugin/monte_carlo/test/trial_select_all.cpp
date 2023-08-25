#include "utils/test/utils.h"
#include "monte_carlo/include/trial_select_all.h"

namespace feasst {

TEST(TrialSelectAll, serialize) {
  TrialSelectAll sel;
  TrialSelectAll sel2 = test_serialize(sel);
}

}  // namespace feasst
