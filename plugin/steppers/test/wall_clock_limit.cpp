#include "utils/test/utils.h"
#include "system/include/utils.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_translate.h"
#include "steppers/include/wall_clock_limit.h"
#include "steppers/include/utils.h"

namespace feasst {

TEST(WallClockLimit, limit) {
  MonteCarlo mc;
  mc.set(lennard_jones());
  mc.set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  mc.add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  add_common_steppers(&mc, {{"steps_per", str(1e4)},
                            {"file_append", "tmp/lj"}});
  TRY(
    mc.add(MakeWallClockLimit({{"max_hours", "1e-9"}}));
    mc.attempt(1e4);
    CATCH_PHRASE("exceed the maximum");
  );
  auto limit = test_serialize(*MakeWallClockLimit({{"max_hours", "1e-9"}}));
}

}  // namespace feasst
