#include <sstream>
#include "utils/test/utils.h"
#include "system/include/thermo_params.h"

namespace feasst {

TEST(ThermoParams, serialize) {
  TRY(
    ThermoParams().beta();
    CATCH_PHRASE("beta must be initialized before use");
  );

  auto params = MakeThermoParams({{"beta", "1.2"}});
  ThermoParams params2 = test_serialize(*params);
  EXPECT_NEAR(params2.beta(), 1.2, NEAR_ZERO);

  TRY(
    params2.pressure();
    CATCH_PHRASE("pressure must be initialized before use");
  );

  auto params3 = MakeThermoParams({{"pressure", "1.2"}});
  ThermoParams params4 = test_serialize(*params3);
  EXPECT_NEAR(params4.pressure(), 1.2, NEAR_ZERO);
}

}  // namespace feasst
