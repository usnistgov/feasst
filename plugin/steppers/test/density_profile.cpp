#include "utils/test/utils.h"
#include "steppers/include/density_profile.h"

namespace feasst {

TEST(DensityProfile, serialize) {
  auto density_profile = MakeDensityProfile({{"output_file", "tmp"}});
  auto density_profile2 = test_serialize<DensityProfile, Analyze>(*density_profile);
}

}  // namespace feasst
