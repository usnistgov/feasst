#include "utils/test/utils.h"
#include "system/include/lennard_jones.h"
#include "system/include/long_range_corrections.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_translate.h"
#include "steppers/include/wall_clock_limit.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/tune.h"
#include "steppers/include/log_and_movie.h"
#include "steppers/include/seek_modify.h"

namespace feasst {

TEST(SeekModify, seek) {
  MonteCarlo mc;
  mc.add(MakeConfiguration({{"cubic_side_length", "8"},
                            {"particle_type0", "../forcefield/lj.fstprt"}}));
  mc.add(MakePotential(MakeLennardJones()));
  mc.add(MakePotential(MakeLongRangeCorrections()));
  mc.set(MakeThermoParams({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc.set(MakeMetropolis());
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  mc.add(MakeLogAndMovie({{"trials_per_write", str(1e4)}, {"file_name", "tmp/lj"}}));
  mc.add(MakeCheckEnergy({{"trials_per_update", str(1e4)}, {"tolerance", str(1e-9)}}));
  mc.add(MakeTune());
  mc.add(MakeTune());
  mc.add(MakeWallClockLimit({{"max_hours", "1e-9"}}));
  EXPECT_EQ(0, SeekModify().index("CheckEnergy", mc)[0]);
  EXPECT_EQ(-1, SeekModify().index("CheckEnergy", mc)[1]);
  EXPECT_EQ(1, SeekModify().index("Tune", mc)[0]);
  EXPECT_EQ(-1, SeekModify().index("Tune", mc)[1]);
  EXPECT_EQ(-1, SeekModify().index("MagicalUnicorn", mc)[0]);
  EXPECT_EQ(-1, SeekModify().index("MagicalUnicorn", mc)[0]);
}

}  // namespace feasst
