#include "utils/test/utils.h"
#include "models/include/model_square_well.h"

namespace feasst {

TEST(ModelSquareWell, serialize) {
  ModelSquareWell model;
  std::shared_ptr<Model> model2 = test_serialize<ModelSquareWell, Model>(model,
    "ModelSquareWell 1 ");
}

}  // namespace feasst
