#include "utils/test/utils.h"
#include "steppers/include/mean_squared_displacement.h"
#include "monte_carlo/test/monte_carlo_test.h"

namespace feasst {

TEST(MeanSquaredDisplacement, msd) {
  MonteCarlo mc;
  mc_lj(&mc);
  mc.seek_num_particles(50);
  mc.add(MakeMeanSquaredDisplacement({
    {"steps_per_update", "10"},
    {"steps_per_write", "100"},
    {"updates_per_origin", "10"},
    {"file_name", "tmp/msd.txt"},
  }));
  mc.attempt(1e3);
  MonteCarlo mc2 = test_serialize(mc);
}

}  // namespace feasst
