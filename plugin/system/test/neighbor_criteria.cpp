#include "utils/test/utils.h"
#include "system/include/neighbor_criteria.h"

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

}  // namespace feasst
