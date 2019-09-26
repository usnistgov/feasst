#include "utils/test/utils.h"
#include "monte_carlo/include/trial_compute_move.h"

namespace feasst {

TEST(TrialComputeMove, serialize) {
  TrialComputeMove add;
//  std::stringstream ss;
//  add.serialize(ss);
//  INFO(ss.str());
//  TrialComputeMove add2(ss);
//  add2.serialize(ss);
//  INFO(ss.str());
  TrialComputeMove add3 = test_serialize(add);
}

}  // namespace feasst
