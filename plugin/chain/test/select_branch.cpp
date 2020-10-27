#include "utils/test/utils.h"
#include "chain/include/select_branch.h"

namespace feasst {

TEST(SelectBranch, serialize) {
  SelectBranch add({{"mobile_site", "2"}, {"mobile_site2", "3"}, {"anchor_site", "1"}, {"anchor_site2", "0"}});
  SelectBranch add3 = test_serialize(add);
}

}  // namespace feasst
