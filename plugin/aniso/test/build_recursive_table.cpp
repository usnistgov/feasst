#include "utils/test/utils.h"
#include "aniso/include/build_recursive_table.h"

namespace feasst {

TEST(BuildRecursiveTable, serialize) {
  BuildRecursiveTable obj;
  BuildRecursiveTable obj2 = test_serialize(obj);
}

}  // namespace feasst
