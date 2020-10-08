#include "utils/test/utils.h"
#include "configuration/include/physical_constants.h"

namespace feasst {

TEST(PhysicalConstants, CODATA2018) {
  CODATA2018 con;
  CODATA2018 con2 = con;
  CODATA2018 con3 = test_serialize(con2);
  EXPECT_NEAR(con3.charge_conversion(), 1389.35457644382, 1e-11);
}

}  // namespace feasst
