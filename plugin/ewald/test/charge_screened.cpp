#include <vector>
#include <tuple>
#include "utils/test/utils.h"
#include "ewald/include/charge_screened.h"
#include "ewald/test/system_example.h"

namespace feasst {

TEST(ChargeScreened, SRSW_refconfig) {
  test_cases(
    { std::make_tuple(MakeCODATA2010(), -4646.8607665992467),
      std::make_tuple(MakeCODATA2014(), -4646.8607561619547),
      std::make_tuple(MakeCODATA2018(), -4646.8607600872092) },
    MakeChargeScreened()
  );
}

}  // namespace feasst
