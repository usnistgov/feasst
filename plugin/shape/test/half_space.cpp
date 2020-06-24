#include "utils/test/utils.h"
#include "utils/include/debug.h"
#include "math/include/random_mt19937.h"
#include "shape/include/half_space.h"

namespace feasst {

TEST(HalfSpace, serialize) {
  HalfSpace half_space({
    {"dimension", "2"},
    {"intersection", "1."},
    {"direction", "1"}});
  TRY(
    HalfSpace half_space({
      {"dimension", "2"},
      {"intersection", "1."},
      {"direction", "0"},
    });
    CATCH_PHRASE("invalid direction");
  );
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

double hamaker_half_plane(const double epsilon,
    const double alpha,
    const double distance) {
  return 2*epsilon*PI/(alpha - 3.)/(alpha - 2.)/std::pow(distance, alpha - 3);
}

TEST(HalfSpace, integrate) {
  const double distance = 1.;
  const double alpha0 = 6.;
  const double epsilon0 = -1.;
  HalfSpace half_space({
    {"dimension", "2"},
    {"intersection", "1."},
    {"direction", "-1"}});
  auto random = MakeRandomMT19937();
  const double inte = half_space.integrate(
    Position({0., 0., 0.}), random.get(), {
      {"alpha0", str(alpha0)},
      {"epsilon0", str(epsilon0)},
      {"max_radius", "10"},
      {"num_shells", "1000"},
      {"points_per_shell", "100"}});
  EXPECT_NEAR(inte, hamaker_half_plane(epsilon0, alpha0, distance), 0.04);
}

}  // namespace feasst
