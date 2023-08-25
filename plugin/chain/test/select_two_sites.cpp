#include "utils/test/utils.h"
#include "chain/include/select_two_sites.h"

namespace feasst {

TEST(SelectTwoSites, serialize) {
  SelectTwoSites add({{"mobile_site", "1"}, {"mobile_site2", "0"}});
  SelectTwoSites add2 = test_serialize(add);
}

}  // namespace feasst
