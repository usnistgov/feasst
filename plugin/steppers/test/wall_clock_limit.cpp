#include "utils/test/utils.h"
#include "system/include/lennard_jones.h"
#include "system/include/long_range_corrections.h"
#include "system/include/potential.h"
#include "system/include/thermo_params.h"
#include "monte_carlo/test/monte_carlo_utils.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_translate.h"
#include "steppers/include/wall_clock_limit.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/tune.h"
#include "steppers/include/log.h"
#include "steppers/include/movie.h"
#include "steppers/include/seek_analyze.h"

namespace feasst {

TEST(WallClockLimit, limit) {
  MonteCarlo mc;
  mc.add(MakeConfiguration({{"cubic_side_length", "8"},
                            {"particle_type0", "../particle/lj.fstprt"}}));
  mc.add(MakePotential(MakeLennardJones()));
  mc.add(MakePotential(MakeLongRangeCorrections()));
  mc.set(MakeThermoParams({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc.add(MakeLog({{"trials_per_write", str(1e4)}, {"output_file", "tmp/lj.txt"}}));
  mc.add(MakeMovie({{"trials_per_write", str(1e4)}, {"output_file", "tmp/lj.xyz"}}));
  mc.add(MakeCheckEnergy({{"trials_per_update", str(1e4)}, {"tolerance", str(1e-9)}}));
  mc.add(MakeTune());
  mc.add(MakeWallClockLimit({{"max_hours", "1e-9"}}));
  TRY(
    mc.attempt(1e4);
    CATCH_PHRASE("exceed the maximum");
  );
  auto limit = test_serialize(*MakeWallClockLimit({{"max_hours", "1e-9"}}));
  EXPECT_EQ(0, SeekAnalyze().index("Log", mc)[0]);
  EXPECT_EQ(-1, SeekAnalyze().index("Log", mc)[1]);
  EXPECT_EQ(1, SeekAnalyze().index("Movie", mc)[0]);
  EXPECT_EQ(-1, SeekAnalyze().index("Movie", mc)[1]);
  EXPECT_EQ(2, SeekAnalyze().index("WallClockLimit", mc)[0]);
  EXPECT_EQ(-1, SeekAnalyze().index("WallClockLimit", mc)[1]);
  EXPECT_EQ(-1, SeekAnalyze().index("MagicalUnicorn", mc)[0]);
  EXPECT_EQ(-1, SeekAnalyze().index("MagicalUnicorn", mc)[0]);
}

}  // namespace feasst
