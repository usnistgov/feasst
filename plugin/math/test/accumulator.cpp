#include "utils/test/utils.h"
#include "math/include/accumulator.h"
#include "math/include/constants.h"

namespace feasst {

TEST(Accumulator, constructor) {
  Accumulator a;
  a.set_block(5);
  EXPECT_NEAR(-NEAR_INFINITY, a.max(), 1);
  EXPECT_NEAR(NEAR_INFINITY, a.min(), 1);

  TRY(
    a.last_value();
    CATCH_PHRASE("no values accumul");
  );
  for (int i = 0; i < 20; ++i) {
    a.accumulate(i);
  }
  EXPECT_EQ(19, a.last_value());
  EXPECT_EQ(20, a.num_values());
  EXPECT_NEAR(19/2., a.average(), NEAR_ZERO);
  EXPECT_NEAR(5.916079783, a.std(), 1e-10);
  EXPECT_NEAR(3.2274861218, a.block_stdev(), 1e-10);
  EXPECT_NEAR(19, a.max(), NEAR_ZERO);
  EXPECT_NEAR(0, a.min(), NEAR_ZERO);

  // copy a using copy constructor and constructor and test
  //Accumulator b(a.num_values(), a.sum(), a.sum_of_squared());

  // serialize
  Accumulator b = test_serialize(a);

  Accumulator c = a;
  EXPECT_EQ(b.num_values(), c.num_values());
  EXPECT_EQ(b.sum(), c.sum());
  EXPECT_EQ(b.sum_of_squared(), c.sum_of_squared());
  EXPECT_EQ(b.average(), c.average());
  EXPECT_EQ(b.std(), c.std());

  // reset a and test
  a.reset();
  EXPECT_NEAR(0., a.average(), NEAR_ZERO);
  EXPECT_EQ(0, a.num_values());
  EXPECT_EQ(0, a.sum());
  EXPECT_EQ(0, a.sum_of_squared());
}

}  // namespace feasst
