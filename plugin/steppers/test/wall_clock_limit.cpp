#include "utils/test/utils.h"
#include "monte_carlo/include/utils.h"
#include "steppers/include/wall_clock_limit.h"
#include "steppers/include/utils.h"

namespace feasst {

TEST(WallClockLimit, limit) {
  MonteCarlo mc;
  lennard_jones(&mc);
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
