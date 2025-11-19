#include "utils/test/utils.h"
#include "confinement/include/model_table_spherical.h"

namespace feasst {

TEST(ModelTableSphere1D, serialize) {
  auto obj = std::make_unique<ModelTableSphere1D>();
  auto obj2 = test_serialize_unique(*obj);
}

}  // namespace feasst
