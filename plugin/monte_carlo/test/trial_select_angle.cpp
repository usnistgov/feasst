#include "utils/test/utils.h"
#include "monte_carlo/include/trial_select_angle.h"

namespace feasst {

TEST(TrialSelectAngle, serialize) {
  TrialSelectAngle add({{"mobile_site", "2"}, {"anchor_site", "1"}, {"anchor_site2", "0"}});
//  std::stringstream ss;
//  add.serialize(ss);
//  INFO(ss.str());
//  TrialSelectAngle add2(ss);
//  add2.serialize(ss);
//  INFO(ss.str());
  TrialSelectAngle add3 = test_serialize(add);
}

}  // namespace feasst
