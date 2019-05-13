#include "utils/test/utils.h"
#include "confinement/include/half_space.h"
#include "utils/include/debug.h"

namespace feasst {

TEST(Shape, HalfSpace) {
  HalfSpace half_space({
    {"dimension", "2"},
    {"intersection", "1."},
    {"direction", "1"},
  });
  try {
    HalfSpace half_space({
      {"dimension", "2"},
      {"intersection", "1."},
      {"direction", "0"},
    });
    CATCH_PHRASE("invalid direction");
  }
  Position point;
  point.set_vector({15, -56.54, 2.});
  EXPECT_TRUE(half_space.is_inside(point));
  point.set_vector({15, -56.54, 1.000000000000001});
  EXPECT_TRUE(half_space.is_inside(point));
  point.set_vector({15, -56.54, 0.9999999});
  EXPECT_FALSE(half_space.is_inside(point));
  point.set_vector({15, -56.54, 1.500000000001});
  EXPECT_TRUE(half_space.is_inside(point, 1.));
  point.set_vector({15, -56.54, 1.4999999999});
  EXPECT_FALSE(half_space.is_inside(point, 1.));

  // serialize
  auto half_space3 = test_serialize<HalfSpace, Shape>(half_space);
  point.set_vector({15, -56.54, 1.4999999999});
  EXPECT_FALSE(half_space3->is_inside(point, 1.));

  auto half_space2 = HalfSpace({
    {"dimension", "2"},
    {"intersection", "3."},
    {"direction", "-1"},
  });
  EXPECT_TRUE(half_space2.is_inside(point));
  point.set_vector({15, -56.54, 3.000000000000001});
  EXPECT_FALSE(half_space2.is_inside(point));
  point.set_vector({15, -56.54, 2.9999999999});
  EXPECT_TRUE(half_space2.is_inside(point));
}

}  // namespace feasst
