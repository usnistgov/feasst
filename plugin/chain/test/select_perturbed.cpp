#include "utils/test/utils.h"
#include "chain/include/select_perturbed.h"

namespace feasst {

TEST(SelectPerturbed, serialize) {
  SelectPerturbed add;
//  std::stringstream ss;
//  add.serialize(ss);
//  INFO(ss.str());
//  SelectPerturbed add2(ss);
//  add2.serialize(ss);
//  INFO(ss.str());
  SelectPerturbed add2 = test_serialize(add);
}

}  // namespace feasst
