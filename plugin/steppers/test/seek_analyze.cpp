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

TEST(SeekAnalyze, seek) {
  auto mc = MakeMonteCarlo({{
    {"Configuration", {{"cubic_side_length", "8"},
                       {"particle_type0", "../particle/lj.txt"}}},
    {"Potential", {{"Model", "LennardJones"}}},
    {"Potential", {{"VisitModel", "LongRangeCorrections"}}},
    {"ThermoParams", {{"beta", "1.2"}, {"chemical_potential", "1."}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"weight", "1."}, {"tunable_param", "1."}}},
    {"Log", {{"trials_per_write", str(1e4)}, {"output_file", "tmp/lj.txt"}}},
    {"Movie", {{"trials_per_write", str(1e4)}, {"output_file", "tmp/lj.xyz"}}},
    {"CheckEnergy", {{"trials_per_update", str(1e4)}, {"tolerance", str(1e-9)}}},
    {"Tune", {{}}},
    {"WallClockLimit", {{"max_hours", "1e-9"}}},
  }}, true);
  EXPECT_EQ(0, SeekAnalyze().index("Log", *mc)[0]);
  EXPECT_EQ(-1, SeekAnalyze().index("Log", *mc)[1]);
  EXPECT_EQ(1, SeekAnalyze().index("Movie", *mc)[0]);
  EXPECT_EQ(-1, SeekAnalyze().index("Movie", *mc)[1]);
  EXPECT_EQ(2, SeekAnalyze().index("WallClockLimit", *mc)[0]);
  EXPECT_EQ(-1, SeekAnalyze().index("WallClockLimit", *mc)[1]);
  EXPECT_EQ(-1, SeekAnalyze().index("MagicalUnicorn", *mc)[0]);
  EXPECT_EQ(-1, SeekAnalyze().index("MagicalUnicorn", *mc)[0]);
}

}  // namespace feasst
