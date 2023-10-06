#include "utils/test/utils.h"
#include "monte_carlo/include/always_reject.h"

namespace feasst {

TEST(AlwaysReject, serialize) {
  AlwaysReject reject;
  AlwaysReject reject2 = test_serialize(reject);
}

}  // namespace feasst
