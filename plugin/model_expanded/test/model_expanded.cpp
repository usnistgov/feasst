#include "utils/test/utils.h"
#include "model_expanded/include/model_expanded.h"

namespace feasst {

TEST(ModelExpanded, serialize) {
  ModelExpanded model;
  std::shared_ptr<Model> model2 = test_serialize<ModelExpanded, Model>(model,
    "ModelExpanded 2094 -1 -1 -1 -1 573 0 7168 0 ");
}

}  // namespace feasst
