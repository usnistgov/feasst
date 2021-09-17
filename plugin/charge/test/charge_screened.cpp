#include <vector>
#include <tuple>
#include "utils/test/utils.h"
#include "charge/include/charge_screened.h"
#include "charge/test/system_example.h"

namespace feasst {

TEST(ChargeScreened, SRSW_refconfig) {
  test_cases(
    { std::make_tuple(MakeCODATA2010(), -4646.8607665992467),
      std::make_tuple(MakeCODATA2014(), -4646.8607561619547),
      std::make_tuple(MakeCODATA2018(), -4646.8607600872092) },
    MakeChargeScreened({{"table_size", "0"}})
  );
  test_cases(
    { std::make_tuple(MakeCODATA2010(), -4646.8608952065752),
      std::make_tuple(MakeCODATA2014(), -4646.8608847693049),
      std::make_tuple(MakeCODATA2018(), -4646.8608886944894)},
    MakeChargeScreened({{"table_size", str(int(1e5))}})
  );
}

}  // namespace feasst
