#include "utils/test/utils.h"
#include "configuration/include/configuration.h"
#include "models/include/two_body_table.h"
#include "aniso/include/model_recursive_table.h"

namespace feasst {

TEST(ModelRecursiveTable, serialize) {
  auto obj = std::make_shared<ModelRecursiveTable>();
  //auto obj = std::make_unique<ModelRecursiveTable>();
  auto obj2 = test_serialize(*obj);
  //auto obj2 = test_serialize_unique(*obj);
}

}  // namespace feasst
