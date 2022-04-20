#include "utils/test/utils.h"
#include "system/include/lennard_jones.h"
#include "system/include/long_range_corrections.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_translate.h"
#include "steppers/include/log.h"
#include "steppers/include/movie.h"

namespace feasst {

TEST(Steppers, CheckForSameFileName) {
  MonteCarlo mc;
  mc.add(MakeConfiguration({{"cubic_box_length", "8"},
                            {"particle_type0", "../forcefield/lj.fstprt"}}));
  mc.add(MakePotential(MakeLennardJones()));
  mc.add(MakePotential(MakeLongRangeCorrections()));
  mc.set(MakeThermoParams({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc.add(MakeMovie({{"trials_per", str(1e4)}, {"file_name", "tmp/lj.xyz"}}));
  try {
    mc.add(MakeLog({{"trials_per", str(1e4)}, {"file_name", "tmp/lj.xyz"}}));
    CATCH_PHRASE("should not have the same file_name");
  }
}

}  // namespace feasst
