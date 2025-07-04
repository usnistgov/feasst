#include "utils/test/utils.h"
#include "system/include/lennard_jones.h"
#include "system/include/long_range_corrections.h"
#include "system/include/potential.h"
#include "system/include/thermo_params.h"
#include "monte_carlo/test/monte_carlo_utils.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_translate.h"
#include "steppers/include/log.h"
#include "steppers/include/movie.h"

namespace feasst {

TEST(Steppers, CheckForSameFileName) {
  MonteCarlo mc;
  mc.add(MakeConfiguration({{"cubic_side_length", "8"},
                            {"particle_type0", "../particle/lj.txt"}}));
  mc.add(MakePotential(MakeLennardJones()));
  mc.add(MakePotential(MakeLongRangeCorrections()));
  mc.set(MakeThermoParams({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc.add(MakeMovie({{"trials_per_write", str(1e4)}, {"output_file", "tmp/lj.xyz"}}));
  TRY(
    mc.add(MakeLog({{"trials_per_write", str(1e4)}, {"output_file", "tmp/lj.xyz"}}));
    CATCH_PHRASE("should not have the same output_file");
  );
}

}  // namespace feasst
