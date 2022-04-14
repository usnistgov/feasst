#include "utils/test/utils.h"
#include "chain/include/select_particle_pivot.h"

namespace feasst {

TEST(SelectParticlePivot, serialize) {
  std::shared_ptr<SelectParticlePivot> sel;
  sel = MakeSelectParticlePivot({{"particle_type", "0"}});
  SelectParticlePivot sel2 = test_serialize(*sel);
}

}  // namespace feasst
