#include "utils/test/utils.h"
#include "configuration/include/model_params.h"

namespace feasst {

TEST(Epsilon, serialize) {
  Epsilon param;
  Epsilon param2 = test_serialize(param);
  auto param3 = std::make_shared<Epsilon>();
  Epsilon param4 = test_serialize(*param3);
}

}
