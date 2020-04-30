#include "utils/test/utils.h"
#include "chain/include/select_segment.h"

namespace feasst {

TEST(SelectSegment, serialize) {
  SelectSegment add;
//  std::stringstream ss;
//  add.serialize(ss);
//  INFO(ss.str());
//  SelectSegment add2(ss);
//  add2.serialize(ss);
//  INFO(ss.str());
  SelectSegment add2 = test_serialize(add);
}

}  // namespace feasst
