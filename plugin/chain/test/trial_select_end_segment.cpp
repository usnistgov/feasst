#include "utils/test/utils.h"
#include "chain/include/trial_select_end_segment.h"

namespace feasst {

TEST(TrialSelectEndSegment, serialize) {
  TrialSelectEndSegment add;
//  std::stringstream ss;
//  add.serialize(ss);
//  INFO(ss.str());
//  TrialSelectEndSegment add2(ss);
//  add2.serialize(ss);
//  INFO(ss.str());
  TrialSelectEndSegment add2 = test_serialize(add);
}

}  // namespace feasst
