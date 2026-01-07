#include "utils/test/utils.h"
#include "configuration/include/configuration.h"
#include "models/include/two_body_table.h"
#include "models/include/recursive_table_potential.h"

namespace feasst {

TEST(RecursiveTablePotential, serialize) {
  auto obj = std::make_unique<RecursiveTablePotential>();
  auto obj2 = test_serialize_unique(*obj);
}

}  // namespace feasst
