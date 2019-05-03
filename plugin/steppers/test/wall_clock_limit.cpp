#include <memory>
#include <gtest/gtest.h>
#include "core/test/monte_carlo_test.h"
#include "steppers/include/wall_clock_limit.h"

namespace feasst {

TEST(WallClockLimit, limit) {
  MonteCarlo mc = mc_lj();
  try {
    mc.add(MakeWallClockLimit({{"max_hours", "1e-9"}}));
    mc.attempt(1e4);
    CATCH_PHRASE("exceed the maximum");
  }
}

}  // namespace feasst
