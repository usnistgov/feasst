#include "utils/test/utils.h"
#include "configuration/include/configuration.h"
#include "system/include/lennard_jones.h"
#include "system/include/long_range_corrections.h"
#include "system/include/potential.h"
#include "system/include/thermo_params.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_translate.h"
#include "steppers/include/wall_clock_limit.h"
#include "steppers/include/check_energy.h"
#include "steppers/include/tune.h"
#include "steppers/include/seek_modify.h"

namespace feasst {

TEST(SeekModify, seek) {
  auto mc = MakeMonteCarlo({{
    {"Configuration", {{"cubic_side_length", "8"},
                       {"particle_type0", "../particle/lj.fstprt"}}},
    {"Potential", {{"Model", "LennardJones"}}},
    {"Potential", {{"VisitModel", "LongRangeCorrections"}}},
    {"ThermoParams", {{"beta", "1.2"}, {"chemical_potential", "1."}}},
    {"Metropolis", {{}}},
    {"TrialTranslate", {{"weight", "1."}, {"tunable_param", "1."}}},
  //  {"LogAndMovie({{"trials_per_write", str(1e4)}, {"output_file", "tmp/lj"}}},
    {"CheckEnergy", {{"trials_per_update", str(1e4)}, {"tolerance", str(1e-9)}}},
    {"Tune", {{}}},
    {"Tune", {{}}},
    {"WallClockLimit", {{"max_hours", "1e-9"}}},
  }}, true);
  EXPECT_EQ(0, SeekModify().index("CheckEnergy", *mc)[0]);
  EXPECT_EQ(-1, SeekModify().index("CheckEnergy", *mc)[1]);
  EXPECT_EQ(1, SeekModify().index("Tune", *mc)[0]);
  EXPECT_EQ(-1, SeekModify().index("Tune", *mc)[1]);
  EXPECT_EQ(-1, SeekModify().index("MagicalUnicorn", *mc)[0]);
  EXPECT_EQ(-1, SeekModify().index("MagicalUnicorn", *mc)[0]);
}

}  // namespace feasst
