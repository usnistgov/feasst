#include "utils/test/utils.h"
#include "beta_expanded/include/select_nothing.h"

namespace feasst {

TEST(SelectNothing, serialize) {
  SelectNothing select;
  SelectNothing select2 = test_serialize(select);
}

}  // namespace feasst
