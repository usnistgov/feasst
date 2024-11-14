#include "utils/test/utils.h"
#include "confinement/include/zero_background.h"

namespace feasst {

TEST(ZeroBackground, serialize) {
  auto obj = MakeZeroBackground();
  ZeroBackground obj2 = test_serialize(*obj);
}

}  // namespace feasst
