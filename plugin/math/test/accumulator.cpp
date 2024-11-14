#include "utils/test/utils.h"
#include "math/include/formula.h"
#include "math/include/accumulator.h"
#include "math/include/constants.h"

namespace feasst {

TEST(Accumulator, constructor) {
  Accumulator a;
  EXPECT_NEAR(-NEAR_INFINITY, a.max(), 1);
  EXPECT_NEAR(NEAR_INFINITY, a.min(), 1);

  TRY(
    a.last_value();
    CATCH_PHRASE("no values accumul");
  );
  std::vector<double> vals;
  for (int i = 0; i < 20; ++i) {
    a.accumulate(i);
    vals.push_back(i);
  }
  EXPECT_EQ(19, a.last_value());
  EXPECT_EQ(20, a.num_values());
  EXPECT_NEAR(19/2., a.average(), NEAR_ZERO);
  EXPECT_NEAR(5.916079783, a.std(), 1e-10);
  EXPECT_NEAR(1.3228756555322954, a.block_stdev(), 1e-10);
  EXPECT_NEAR(19, a.max(), NEAR_ZERO);
  EXPECT_NEAR(0, a.min(), NEAR_ZERO);

//  INFO("block size " << a.block_size()[0]);
//  for (const auto bl : a.blocks()) {
//    INFO(feasst_str(bl));
//  }
  EXPECT_EQ(vals, a.blocks()[0]);

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

TEST(Accumulator, is_equivalent) {
  auto a = MakeAccumulator();
  auto b = MakeAccumulator();
  std::vector<double> avals = {5, 6, 4, 5.5, 4, 6, 5.1};
  std::vector<double> bvals = {8, 9, 7, 6.8, 10, 5};
  for (const double v : avals) a->accumulate(v);
  for (const double v : bvals) b->accumulate(v);
  EXPECT_FALSE(a->is_equivalent(*b, 2, 1));
  EXPECT_TRUE(a->is_equivalent(*b, 10, 1));
}

class Square : public Formula {
 public:
  double evaluate(const double y) const override { return y*y; }
};

TEST(Accumulator, largest_block) {
  auto a = MakeAccumulator();
  std::vector<double> avals = {5, 6, 4, 5.5, 4, 6, 5.1};
  for (const double v : avals) a->accumulate(v);
  EXPECT_EQ(a->blocks()[0].size(), 7);
  EXPECT_EQ(a->blocks()[1].size(), 3);
  EXPECT_EQ(a->blocks()[2].size(), 1);
  EXPECT_EQ(a->blocks()[3].size(), 0);
  EXPECT_EQ(a->blocks()[0][0], 5);
  EXPECT_EQ(a->blocks()[0][1], 6);
  EXPECT_EQ(a->blocks()[0][2], 4);
  EXPECT_EQ(a->blocks()[0][3], 5.5);
  EXPECT_EQ(a->blocks()[1][0], 5.5);
  EXPECT_EQ(a->blocks()[1][1], 4.75);
  EXPECT_EQ(a->blocks()[1][2], 5);
  EXPECT_EQ(a->blocks()[2][0], 5.125);
  EXPECT_EQ(a->largest_block_operation(), 1);
  EXPECT_EQ(a->largest_blocks()[0], 5.5);
  EXPECT_EQ(a->largest_blocks()[1], 4.75);
  EXPECT_EQ(a->largest_blocks()[2], 5);
  EXPECT_EQ(a->largest_blocks().size(), 3);
  Square sq;
  INFO(a->block_stdev(sq));
}

TEST(Accumulator, serialize_with_inf) {
  Accumulator acc;
  acc.accumulate(std::numeric_limits<long double>::max());
  acc.accumulate(std::numeric_limits<long double>::max());
//  Accumulator acc2 = test_serialize(acc);
}

}  // namespace feasst
