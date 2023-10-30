#include <fstream>
#include "utils/test/utils.h"
#include "utils/include/checkpoint.h"

namespace feasst {

// If this test fails, it may be optimizing the product calculation too quickly
// to checkpoint.
TEST(Checkpoint, init) {
  Checkpoint check({{"checkpoint_file", "tmp/checkpoint"}, {"num_hours", "1e-7"} });
  double prod = 1;
  for (int i = 0; i < 1000; ++i) {
    prod *= i;
    check.check(check);
  }
  std::stringstream ss;
  check.serialize(ss);
  Checkpoint check2(ss);
  EXPECT_EQ(check2.num_hours(), 1e-7);
  test_serialize(check);

  check2.write(check2);
  Checkpoint check3;
  MakeCheckpoint({{"checkpoint_file", "tmp/checkpoint"}})->read(&check3);
  EXPECT_EQ(check3.num_hours(), 1e-7);
}

}  // namespace feasst
