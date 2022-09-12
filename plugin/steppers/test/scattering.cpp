#include "utils/test/utils.h"
#include "steppers/include/scattering.h"

namespace feasst {

TEST(Scattering, serialize) {
  auto an = MakeScattering();
  auto an2 = test_serialize<Scattering, Analyze>(*an);
}

}  // namespace feasst
