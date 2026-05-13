#include "utils/test/utils.h"
#include "steppers/include/extensive_moments.h"

namespace feasst {

TEST(ExtensiveMoments, serialize) {
  auto an = MakeExtensiveMoments();
  auto an2 = test_serialize<ExtensiveMoments, Analyze>(*an);
}

}  // namespace feasst
