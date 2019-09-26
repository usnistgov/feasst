#include "utils/test/utils.h"
#include "chain/include/trial_select_reptate.h"

namespace feasst {

TEST(TrialSelectReptate, serialize) {
  TrialSelectReptate add;
//  std::stringstream ss;
//  add.serialize(ss);
//  INFO(ss.str());
//  TrialSelectReptate add2(ss);
//  add2.serialize(ss);
//  INFO(ss.str());
  TrialSelectReptate add2 = test_serialize(add);
}

}  // namespace feasst
