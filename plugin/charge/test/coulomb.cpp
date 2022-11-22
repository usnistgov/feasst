#include "utils/test/utils.h"
#include "charge/include/coulomb.h"

namespace feasst {

TEST(Coulomb, serialize) {
  Coulomb model;
  std::shared_ptr<Model> model2 = test_serialize<Coulomb, Model>(model,
    "Coulomb 2094 -1 -1 -1 -1 1634 0 ");
}

}  // namespace feasst
