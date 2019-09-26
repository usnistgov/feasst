#include "utils/test/utils.h"
#include "chain/include/trial_select_segment.h"

namespace feasst {

TEST(TrialSelectSegment, serialize) {
  TrialSelectSegment add;
//  std::stringstream ss;
//  add.serialize(ss);
//  INFO(ss.str());
//  TrialSelectSegment add2(ss);
//  add2.serialize(ss);
//  INFO(ss.str());
  TrialSelectSegment add2 = test_serialize(add);
}

}  // namespace feasst
