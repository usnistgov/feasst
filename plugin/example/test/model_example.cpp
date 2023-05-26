#include "utils/test/utils.h"
#include "configuration/include/model_params.h"
#include "example/include/model_example.h"

namespace feasst {

TEST(ModelExample, model_example) {
  auto example = MakeModelExample({{"example_argument", "0.5"}});
  std::shared_ptr<Model> model = test_serialize<ModelExample, Model>(*example);
}

}  // namespace feasst
