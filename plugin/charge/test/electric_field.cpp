#include "utils/test/utils.h"
#include "charge/include/electric_field.h"

namespace feasst {

TEST(ElectricField, serialize) {
  auto model = MakeElectricField({{"field_strength", "1"}});
  std::shared_ptr<Model> model2 = test_serialize<ElectricField, Model>(*model);
}

}  // namespace feasst
