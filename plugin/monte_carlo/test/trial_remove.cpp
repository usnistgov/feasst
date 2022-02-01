#include "utils/test/utils.h"
#include "monte_carlo/include/trial_remove.h"

namespace feasst {

TEST(TrialRemove, serialize) {
  auto remove = MakeTrialRemove();
//  std::stringstream ss;
//  add.serialize(ss);
//  INFO(ss.str());
//  TrialRemove add2(ss);
//  add2.serialize(ss);
//  INFO(ss.str());
  Trial remove2 = test_serialize(*remove);
}

}  // namespace feasst
