#include "utils/test/utils.h"
#include "rxvb/include/select_anchor.h"

namespace feasst {

TEST(SelectAnchor, serialize) {
  SelectAnchor sel;
  SelectAnchor sel2 = test_serialize(sel);
}

}  // namespace feasst
