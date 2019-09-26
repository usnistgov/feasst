#include "utils/test/utils.h"
#include "monte_carlo/include/trial_add.h"

namespace feasst {

TEST(TrialAdd, serialize) {
  TrialAdd add;
//  std::stringstream ss;
//  add.serialize(ss);
//  INFO(ss.str());
//  TrialAdd add2(ss);
//  add2.serialize(ss);
//  INFO(ss.str());
  TrialAdd add2 = test_serialize(add);
}

}  // namespace feasst
