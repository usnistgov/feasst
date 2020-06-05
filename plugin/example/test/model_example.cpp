#include "utils/test/utils.h"
#include "configuration/include/model_params.h"
#include "example/include/model_example.h"

namespace feasst {

TEST(ModelExample, model_example) {
  ModelExample example;

  // test energy
  EXPECT_EQ(0., example.energy(0, 0, 0, ModelParams()));

  std::shared_ptr<Model> model = test_serialize<ModelExample, Model>(example);
}

}  // namespace feasst
