#include "utils/test/utils.h"
#include "monte_carlo/include/trial_compute_add.h"

namespace feasst {

TEST(TrialComputeAdd, serialize) {
  TrialComputeAdd add;
//  std::stringstream ss;
//  add.serialize(ss);
//  INFO(ss.str());
//  TrialComputeAdd add2(ss);
//  add2.serialize(ss);
//  INFO(ss.str());
  TrialComputeAdd add3 = test_serialize(add);
}

}  // namespace feasst
