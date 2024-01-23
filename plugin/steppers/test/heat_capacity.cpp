#include "utils/test/utils.h"
#include "steppers/include/heat_capacity.h"

namespace feasst {

TEST(HeatCapacity, serialize) {
  auto an = MakeHeatCapacity();
  auto an2 = test_serialize<HeatCapacity, Analyze>(*an);
}

}  // namespace feasst
