#include "utils/test/utils.h"
#include "steppers/include/mean_squared_displacement.h"
#include "steppers/include/utils.h"
#include "monte_carlo/include/utils.h"
#include "monte_carlo/include/seek_num_particles.h"

namespace feasst {

TEST(MeanSquaredDisplacement, msd) {
  MonteCarlo mc;
  lennard_jones(&mc);
  add_common_steppers(&mc, {{"steps_per", str(1e4)},
                            {"file_append", "tmp/lj"}});
  SeekNumParticles(50).with_trial_add().run(&mc);
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
