#include <sstream>
#include "utils/test/utils.h"
#include "system/include/model_hard_sphere.h"

namespace feasst {

TEST(ModelHardSphere, serialize) {
  ModelHardSphere model;
  std::shared_ptr<Model> model2 = test_serialize<ModelHardSphere, Model>(model,
    "ModelHardSphere 607 ");
}

}  // namespace feasst
