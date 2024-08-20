#include "utils/test/utils.h"
#include "flat_histogram/include/wang_landau.h"

namespace feasst {

TEST(WangLandau, args) {
  TRY(
    MakeWangLandau();
    CATCH_PHRASE("key(min_flatness) is required");
  );
}

TEST(WangLandau, serialize) {
  auto bias = MakeWangLandau({{"min_flatness", "2"}});
  auto bias2 = test_serialize_unique(*bias);
}

}  // namespace feasst
