#include "utils/test/utils.h"
#include "monte_carlo/test/monte_carlo_test.h"
#include "prefetch/include/prefetch.h"

namespace feasst {

TEST(Prefetch, NVT_benchmark) {
  Prefetch mc;
  mc_lj(&mc);
  mc.set(MakeRandomMT19937({{"seed", "default"}}));
  mc.activate_prefetch(false);
  mc.seek_num_particles(50);
  // //activate prefetch after initial configuration
  mc.activate_prefetch(true);
  // mc.attempt(1e6);  // ~3.5 seconds (now 4.1)
  mc.attempt(1e2);
  //INFO("num " << mc.trials().num_attempts());
}

}  // namespace feasst
