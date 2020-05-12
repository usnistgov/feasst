#include "utils/test/utils.h"
#include "flat_histogram/include/wang_landau.h"

namespace feasst {

TEST(WangLandau, args) {
  TRY(
    MakeWangLandau();
    CATCH_PHRASE("key(min_flatness) is required");
  );
}

}  // namespace feasst
