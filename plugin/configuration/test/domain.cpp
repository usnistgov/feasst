#include "utils/test/utils.h"
#include "configuration/include/domain.h"
#include "math/include/constants.h"
#include "math/include/random_mt19937.h"

namespace feasst {

TEST(Domain, Volume) {
  Domain domain;
  domain.set_cubic(5);
  EXPECT_EQ(125, domain.volume());
}

TEST(Domain, min_side_length) {
  Domain domain;
  domain.set_cubic(10);
  EXPECT_NEAR(10, domain.min_side_length(), NEAR_ZERO);
  domain.set_side_length(1, 5.234);
  EXPECT_NEAR(5.234, domain.min_side_length(), NEAR_ZERO);
}

// compute the volume of a sphere
TEST(Domain, random_position) {
  Domain domain;
  RandomMT19937 random;
  domain.set_cubic(2);
  int inside = 0;
  const int trials = 1e3;
  for (int trial = 0; trial < trials; ++trial) {
    if (domain.random_position(&random).squared_distance() <= 1) {
      ++inside;
    }
  }
  EXPECT_NEAR(static_cast<double>(inside)/trials, PI/6, 0.06);
}

TEST(Domain, wrap) {
  auto domain = MakeDomain({{"cubic_box_length", "5"}});
  Position pos;
  pos.set_vector({5, 5, 5});
  Position shift = domain->shift(pos);
  EXPECT_FALSE(domain->is_tilted());
  EXPECT_NEAR(-5, shift.coord(0), NEAR_ZERO);
  EXPECT_NEAR(-5, shift.coord(1), NEAR_ZERO);
  EXPECT_NEAR(-5, shift.coord(2), NEAR_ZERO);

  domain = MakeDomain({{"cubic_box_length", "5"},
    {"xy", "1"}});
  shift = domain->shift(pos);
  EXPECT_TRUE(domain->is_tilted());
  EXPECT_NEAR(-6, shift.coord(0), NEAR_ZERO);
  EXPECT_NEAR(-5, shift.coord(1), NEAR_ZERO);
  EXPECT_NEAR(-5, shift.coord(2), NEAR_ZERO);

  domain = MakeDomain({{"cubic_box_length", "5"},
    {"xy", "1"}, {"xz", "1"}});
  shift = domain->shift(pos);
  EXPECT_TRUE(domain->is_tilted());
  EXPECT_NEAR(-7, shift.coord(0), NEAR_ZERO);
  EXPECT_NEAR(-5, shift.coord(1), NEAR_ZERO);
  EXPECT_NEAR(-5, shift.coord(2), NEAR_ZERO);

  domain = MakeDomain({{"cubic_box_length", "5"},
    {"xy", "1"}, {"xz", "1"}, {"yz", "1"}});
  shift = domain->shift(pos);
  EXPECT_TRUE(domain->is_tilted());
  EXPECT_NEAR(-7, shift.coord(0), NEAR_ZERO);
  EXPECT_NEAR(-6, shift.coord(1), NEAR_ZERO);
  EXPECT_NEAR(-5, shift.coord(2), NEAR_ZERO);

  // serialize
  Domain domain2 = test_serialize(*domain);
  EXPECT_EQ(domain->volume(), domain2.volume());
  EXPECT_EQ(domain->yz(), domain2.yz());
  EXPECT_EQ(1., domain2.yz());
  EXPECT_EQ(domain->periodic(2), domain2.periodic(2));
}

TEST(Domain, non_cubic) {
  auto domain = MakeDomain({
    {"side_length0", "3"},
    {"side_length1", "4"},
    {"side_length2", "5"},
  });
  EXPECT_EQ(domain->dimension(), 3);
  EXPECT_EQ(domain->volume(), 3*4*5);
  Position shift = domain->shift_opt(Position({4.5, 2, 1}));
}

}  // namespace feasst
