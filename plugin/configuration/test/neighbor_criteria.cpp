#include "utils/test/utils.h"
#include "configuration/include/neighbor_criteria.h"

namespace feasst {

TEST(NeighborCriteria, volume) {
  NeighborCriteria neighcrit({
    {"minimum_distance", "0.8"},
    {"maximum_distance", "3."},
    {"site_type0", "0"},
    {"site_type1", "0"},
  });
  EXPECT_NEAR(neighcrit.volume(3), 110.95267494438191, NEAR_ZERO);
  EXPECT_NEAR(neighcrit.volume(2), 26.263714584010668, NEAR_ZERO);
  test_serialize(neighcrit);
  TRY(
    neighcrit.volume(1);
    CATCH_PHRASE("dimension: 1");
  );
}

TEST(NeighborCriteria, type) {
  NeighborCriteria neighcrit({
    {"minimum_distance", "0.8"},
    {"maximum_distance", "3."},
    {"site_type0", "0"},
    {"site_type1", "1"},
  });
  EXPECT_FALSE(neighcrit.is_accepted(0, 1, 0, 0));
  EXPECT_TRUE(neighcrit.is_accepted(0, 1, 0, 1));
  NeighborCriteria nc2({
    {"minimum_distance", "0.8"},
    {"maximum_distance", "3."},
    {"site_type0", "0"},
    {"site_type1", "1"},
    {"site_type0_alt", "2"},
    {"site_type1_alt", "3"},
  });
  auto nc3 = test_serialize(nc2);
  EXPECT_FALSE(nc3.is_accepted(0, 1, 0, 0));
  EXPECT_TRUE (nc3.is_accepted(0, 1, 0, 1));
  EXPECT_FALSE(nc3.is_accepted(0, 1, 0, 2));
  EXPECT_FALSE(nc3.is_accepted(0, 1, 0, 3));
  EXPECT_TRUE (nc3.is_accepted(0, 1, 1, 0));
  EXPECT_FALSE(nc3.is_accepted(0, 1, 1, 1));
  EXPECT_FALSE(nc3.is_accepted(0, 1, 1, 2));
  EXPECT_FALSE(nc3.is_accepted(0, 1, 1, 3));
  EXPECT_FALSE(nc3.is_accepted(0, 1, 2, 0));
  EXPECT_FALSE(nc3.is_accepted(0, 1, 2, 1));
  EXPECT_FALSE(nc3.is_accepted(0, 1, 2, 2));
  EXPECT_TRUE (nc3.is_accepted(0, 1, 2, 3));
  EXPECT_FALSE(nc3.is_accepted(0, 1, 3, 0));
  EXPECT_FALSE(nc3.is_accepted(0, 1, 3, 1));
  EXPECT_TRUE (nc3.is_accepted(0, 1, 3, 2));
  EXPECT_FALSE(nc3.is_accepted(0, 1, 3, 3));
}

}  // namespace feasst
