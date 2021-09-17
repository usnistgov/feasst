#include "utils/test/utils.h"
#include "charge/include/charge_self.h"
#include "charge/test/system_example.h"

namespace feasst {

TEST(ChargeSelf, SRSW_refconfig) {
  test_cases(
    { std::make_tuple(MakeCODATA2010(), -23652.08040365018),
      std::make_tuple(MakeCODATA2014(), -23652.080350525335),
      std::make_tuple(MakeCODATA2018(), -23652.080370504391) },
    MakeChargeSelf()
  );
}

}  // namespace feasst
