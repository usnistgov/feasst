#include "utils/test/utils.h"
#include "rxvb/include/select_two_particles.h"

namespace feasst {

TEST(SelectTwoParticles, serialize) {
  TRY(
    auto sel = std::make_unique<SelectTwoParticles>(argtype({{"particle_types", "1,1"}}));
    CATCH_PHRASE("appears more than once");
  );
  auto sel = std::make_unique<SelectTwoParticles>(argtype({{"particle_types", "1,2"}}));
  auto sel2 = test_serialize(sel);
}

}  // namespace feasst
