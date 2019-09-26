#include "utils/test/utils.h"
#include "monte_carlo/include/trial_select_bond.h"

namespace feasst {

TEST(TrialSelectBond, serialize) {
  TrialSelectBond add({{"mobile_site", "1"}, {"anchor_site", "0"}});
//  std::stringstream ss;
//  add.serialize(ss);
//  INFO(ss.str());
//  TrialSelectBond add2(ss);
//  add2.serialize(ss);
//  INFO(ss.str());
  TrialSelectBond add2 = test_serialize(add);
}

}  // namespace feasst
