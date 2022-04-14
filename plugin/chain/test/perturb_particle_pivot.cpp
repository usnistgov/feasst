#include "utils/test/utils.h"
#include "chain/include/perturb_particle_pivot.h"

namespace feasst {

TEST(PerturbParticlePivot, serialize) {
  PerturbParticlePivot pert;
  PerturbParticlePivot pert2 = test_serialize(pert);
}

}  // namespace feasst
