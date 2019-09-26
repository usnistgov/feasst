#include "utils/test/utils.h"
#include "chain/include/trial_select_perturbed.h"

namespace feasst {

TEST(TrialSelectPerturbed, serialize) {
  TrialSelectPerturbed add;
//  std::stringstream ss;
//  add.serialize(ss);
//  INFO(ss.str());
//  TrialSelectPerturbed add2(ss);
//  add2.serialize(ss);
//  INFO(ss.str());
  TrialSelectPerturbed add2 = test_serialize(add);
}

}  // namespace feasst
