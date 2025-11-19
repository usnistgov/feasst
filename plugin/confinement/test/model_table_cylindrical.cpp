#include "utils/test/utils.h"
#include "confinement/include/model_table_cylindrical.h"

namespace feasst {

TEST(ModelTableCylinder1D, serialize) {
  auto obj = std::make_unique<ModelTableCylinder1D>(argtype({{"radius", "1"},
  {"first_point", "0,0,0"}, {"second_point", "0,0,0"}}));
  auto obj2 = test_serialize_unique(*obj);
}

}  // namespace feasst
