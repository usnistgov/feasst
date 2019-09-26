#include "utils/test/utils.h"
#include "monte_carlo/include/trial_compute_remove.h"

namespace feasst {

TEST(TrialComputeRemove, serialize) {
  TrialComputeRemove add;
//  std::stringstream ss;
//  add.serialize(ss);
//  INFO(ss.str());
//  TrialComputeRemove add2(ss);
//  add2.serialize(ss);
//  INFO(ss.str());
  TrialComputeRemove add3 = test_serialize(add);
}

}  // namespace feasst
