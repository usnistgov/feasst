#include "utils/test/utils.h"
#include "chain/include/trial_select_reptate.h"

namespace feasst {

TEST(TrialSelectReptate, serialize) {
  auto add = MakeTrialSelectReptate({{"max_length", "1"}});
//  std::stringstream ss;
//  add.serialize(ss);
//  INFO(ss.str());
//  TrialSelectReptate add2(ss);
//  add2.serialize(ss);
//  INFO(ss.str());
  TrialSelectReptate add2 = test_serialize(*add);
}

}  // namespace feasst
