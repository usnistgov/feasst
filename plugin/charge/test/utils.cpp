
#include "utils/test/utils.h"
#include "configuration/include/configuration.h"
#include "charge/include/utils.h"
#include "system/include/system.h"
#include "charge/test/charge_utils.h"

namespace feasst {

TEST(Charge, utils) {
  System sys = spce({{"physical_constants", "CODATA2010"}, {"alpha", str(5.6/20)}, {"kmax_squared", "38"}});
  EXPECT_EQ(sys.configuration().model_params().physical_constants().class_name(), "CODATA2010");
}

}  // namespace feasst
