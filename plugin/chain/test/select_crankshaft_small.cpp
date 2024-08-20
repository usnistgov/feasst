#include "utils/test/utils.h"
#include "configuration/include/select.h"
#include "chain/include/select_crankshaft_small.h"

namespace feasst {

TEST(SelectCrankshaftSmall, serialize) {
  std::shared_ptr<SelectCrankshaftSmall> sel;
  TRY(
    sel = MakeSelectCrankshaftSmall();
    CATCH_PHRASE("site argument required");
  );
  TRY(
    sel = MakeSelectCrankshaftSmall({{"site", "0"}});
    CATCH_PHRASE("key(anchor_site0) is required");
  );
  TRY(
    sel = MakeSelectCrankshaftSmall({{"site", "0"}, {"anchor_site0", "1"}});
    CATCH_PHRASE("key(anchor_site1) is required");
  );
  sel = MakeSelectCrankshaftSmall({{"site", "0"}, {"site1", "3"}, {"site2", "4"}, {"anchor_site0", "1"}, {"anchor_site1", "2"}});
  SelectCrankshaftSmall sel2 = test_serialize(*sel);
  EXPECT_EQ(sel2.mobile().num_sites(), 3);
  EXPECT_EQ(sel2.anchor().num_sites(), 2);
}

}  // namespace feasst
