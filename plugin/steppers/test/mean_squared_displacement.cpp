#include "utils/test/utils.h"
#include "system/include/utils.h"
#include "steppers/include/mean_squared_displacement.h"
#include "steppers/include/check_energy_and_tune.h"
#include "steppers/include/log_and_movie.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/run.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trials.h"

namespace feasst {

TEST(MeanSquaredDisplacement, msd) {
  MonteCarlo mc;
  mc.set(lennard_jones());
  mc.set(MakeThermoParams({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc.add(MakeTrialAdd({{"particle_type", "0"}}));
  mc.run(MakeRun({{"until_num_particles", "50"}}));
  mc.run(MakeRemoveTrial({{"name", "TrialAdd"}}));
  mc.add(MakeLogAndMovie({{"steps_per", str(1e4)}, {"file_name", "tmp/lj"}}));
  mc.add(MakeCheckEnergyAndTune({{"steps_per", str(1e4)}, {"tolerance", str(1e-9)}}));
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
