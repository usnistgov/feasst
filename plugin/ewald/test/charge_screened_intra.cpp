#include "utils/test/utils.h"
#include "ewald/include/charge_screened_intra.h"
#include "configuration/test/configuration_test.h"
#include "system/include/visit_model_intra.h"
#include "ewald/test/system_example.h"

namespace feasst {

TEST(ChargeScreenedIntra, SRSW_refconfig) {
  test_cases(
    { std::make_tuple(MakeCODATA2010(), 23363.573774608001),
      std::make_tuple(MakeCODATA2014(), 23363.573722131176),
      std::make_tuple(MakeCODATA2018(), 23363.573741866534) },
    MakeChargeScreenedIntra(),
    MakeVisitModelIntra({{"cutoff", "0"}})
  );
}

}  // namespace feasst
