#include <fstream>
#include "utils/test/utils.h"
#include "monte_carlo/include/monte_carlo.h"
#include "actions/include/wait.h"

namespace feasst {

TEST(Wait, serialize) {
  auto action = std::make_shared<Wait>();
  Action action2 = test_serialize(*action);
}

TEST(Wait, parse) {
  std::ofstream file("tmp/tmp");
  MonteCarlo mc;
  Wait(argtype({{"until_file_exists", "tmp/tmp"}})).run(&mc);

  std::ofstream file0("tmp/prefix0suffix");
  std::ofstream file1("tmp/prefix1suffix");
  auto wait = std::make_shared<Wait>(argtype({{"until_file_exists", "tmp/prefix[2]suffix"}}));
  Wait wait2 = test_serialize(*wait);
  EXPECT_EQ(wait2.prefix(), "tmp/prefix");
  EXPECT_EQ(wait2.suffix(), "suffix");
  EXPECT_EQ(wait2.number(), 2);
  wait2.run(&mc);
}

}  // namespace feasst
