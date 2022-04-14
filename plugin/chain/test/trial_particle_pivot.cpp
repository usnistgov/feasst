#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "configuration/include/domain.h"
#include "configuration/include/file_xyz.h"
#include "system/include/system.h"
#include "system/include/lennard_jones.h"
#include "chain/include/trial_particle_pivot.h"

namespace feasst {

TEST(TrialParticlePivot, serialize) {
  auto pivot = MakeTrialParticlePivot({{"particle_type", "0"}});
  test_serialize(*pivot);
}

}  // namespace feasst
