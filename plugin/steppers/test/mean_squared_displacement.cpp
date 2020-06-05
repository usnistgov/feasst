#include "utils/test/utils.h"
#include "steppers/include/mean_squared_displacement.h"
#include "steppers/include/utils.h"
#include "system/include/utils.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/seek_num_particles.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_translate.h"

namespace feasst {

TEST(MeanSquaredDisplacement, msd) {
  MonteCarlo mc;
  mc.set(lennard_jones());
  mc.set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
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
