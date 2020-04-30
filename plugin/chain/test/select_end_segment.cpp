#include "utils/test/utils.h"
#include "chain/include/select_end_segment.h"

namespace feasst {

TEST(SelectEndSegment, serialize) {
  SelectEndSegment add;
//  std::stringstream ss;
//  add.serialize(ss);
//  INFO(ss.str());
//  SelectEndSegment add2(ss);
//  add2.serialize(ss);
//  INFO(ss.str());
  SelectEndSegment add2 = test_serialize(add);
}

}  // namespace feasst
