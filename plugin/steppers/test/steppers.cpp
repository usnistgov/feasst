#include "utils/test/utils.h"
#include "system/include/utils.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trials.h"
#include "steppers/include/log.h"
#include "steppers/include/movie.h"

namespace feasst {

TEST(Steppers, CheckForSameFileName) {
  MonteCarlo mc;
  mc.set(lennard_jones());
  mc.set(MakeThermoParams({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc.add(MakeMovie({{"steps_per", str(1e4)}, {"file_name", "tmp/lj.xyz"}}));
  try {
    mc.add(MakeLog({{"steps_per", str(1e4)}, {"file_name", "tmp/lj.xyz"}}));
    CATCH_PHRASE("should not have the same file_name");
  }
}

}  // namespace feasst
