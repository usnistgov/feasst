#include <memory>
#include "utils/test/utils.h"
#include "monte_carlo/include/metropolis.h"

namespace feasst {

TEST(Metropolis, serialize) {
  //auto crit = MakeMetropolis({{"num_attempts_per_cycle", "1"}});
  //auto crit = MakeMetropolis({{"trials_per_cycle", "1"}, {"num_attempts_per_cycle", "1"}});
  auto crit = MakeMetropolis({{"trials_per_cycle", "1"}});
  Metropolis crit2 = test_serialize(*crit);
}

}  // namespace feasst
