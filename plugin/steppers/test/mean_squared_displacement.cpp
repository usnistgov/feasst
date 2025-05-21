#include "utils/test/utils.h"
#include "monte_carlo/test/monte_carlo_utils.h"
#include "math/include/accumulator.h"
#include "configuration/include/configuration.h"
#include "system/include/lennard_jones.h"
#include "system/include/long_range_corrections.h"
#include "system/include/potential.h"
#include "system/include/thermo_params.h"
#include "steppers/include/mean_squared_displacement.h"
#include "steppers/include/tune.h"
#include "steppers/include/check_energy.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/run.h"
#include "monte_carlo/include/remove.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_add.h"
#include "monte_carlo/include/trial_translate.h"

namespace feasst {

TEST(MeanSquaredDisplacement, msd) {
  MonteCarlo mc;
  mc.add(MakeConfiguration({{"cubic_side_length", "8"},
                            {"particle_type0", "../particle/lj.txt"}}));
  mc.add(MakePotential(MakeLennardJones()));
  mc.add(MakePotential(MakeLongRangeCorrections()));
  mc.set(MakeThermoParams({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc.add(MakeTrialAdd({{"particle_type", "0"}}));
  mc.run(MakeRun({{"until_num_particles", "50"}}));
  mc.run(MakeRemove({{"name", "TrialAdd"}}));
//  mc.add(MakeLogAndMovie({{"trials_per_write", str(1e4)}, {"output_file", "tmp/lj"}}));
  mc.add(MakeCheckEnergy({{"trials_per_update", str(1e4)}, {"tolerance", str(1e-9)}}));
  mc.add(MakeTune());
  mc.add(MakeMeanSquaredDisplacement({
    {"trials_per_update", "10"},
    {"trials_per_write", "100"},
    {"updates_per_origin", "10"},
    {"output_file", "tmp/msd.txt"},
  }));
  mc.attempt(1e3);
  auto mc2 = test_serialize_unique(mc);
}

}  // namespace feasst
