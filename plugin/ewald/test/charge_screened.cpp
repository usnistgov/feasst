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
    MakeChargeScreened({{"disable_table", "true"}})
  );
  test_cases(
    { std::make_tuple(MakeCODATA2010(), -4646.8607667040897),
      std::make_tuple(MakeCODATA2014(), -4646.8607562668349),
      std::make_tuple(MakeCODATA2018(), -4646.8607601920548)},
    MakeChargeScreened({{"disable_table", "false"}})
  );
}

}  // namespace feasst
