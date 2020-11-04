#include <sstream>
#include "utils/test/utils.h"
#include "system/include/hard_sphere.h"

namespace feasst {

TEST(HardSphere, serialize) {
  HardSphere model;
  std::shared_ptr<Model> model2 =
    test_serialize<HardSphere, Model>(model, "HardSphere 607 ");
}

}  // namespace feasst
