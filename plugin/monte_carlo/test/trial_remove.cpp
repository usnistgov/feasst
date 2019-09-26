#include "utils/test/utils.h"
#include "monte_carlo/include/trial_remove.h"

namespace feasst {

TEST(TrialRemove, serialize) {
  TrialRemove add;
//  std::stringstream ss;
//  add.serialize(ss);
//  INFO(ss.str());
//  TrialRemove add2(ss);
//  add2.serialize(ss);
//  INFO(ss.str());
  TrialRemove add2 = test_serialize(add);
}

}  // namespace feasst
