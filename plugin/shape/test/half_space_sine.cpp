#include <fstream>
#include "utils/test/utils.h"
#include "utils/include/debug.h"
#include "math/include/random_mt19937.h"
#include "shape/include/half_space_sine.h"

namespace feasst {

TEST(HalfSpaceSine, args) {
  TRY(
    HalfSpaceSine half_space(
      MakeFormulaSineWave({{"shift", "1"}}),
      { {"dimension", "0"},
        {"direction", "1"},
        {"intersection", "0"},
        {"wave_dimension", "1"}});
    CATCH_PHRASE("use the intersection argument instead of shift");
  );
}

TEST(HalfSpaceSine, direction1) {
  Position point;
  HalfSpaceSine half_space(
    MakeFormulaSineWave({{"amplitude", "0"}}),
    { {"dimension", "0"},
      {"intersection", "0"},
      {"direction", "1"},
      {"wave_dimension", "1"}});
  HalfSpaceSine half_space2 = test_serialize(half_space);
  point.set_vector({0.001, 0, 0});
  EXPECT_TRUE(half_space2.is_inside(point));
  point.set_vector({-0.001, 0, 0});
  EXPECT_FALSE(half_space2.is_inside(point));

  // finite size
  const double diameter = 1.;
  point.set_vector({0.449, 0, 0});
  EXPECT_FALSE(half_space2.is_inside(point, diameter));
  point.set_vector({0.501, 0, 0});
  EXPECT_TRUE(half_space2.is_inside(point, diameter));
}

TEST(HalfSpaceSine, direction2) {
  Position point;
  HalfSpaceSine half_space(
    MakeFormulaSineWave({{"amplitude", "0"}}),
    { {"dimension", "0"},
      {"direction", "-1"},
      {"intersection", "0"},
      {"wave_dimension", "1"}});
  point.set_vector({0.001, 0, 0});
  EXPECT_FALSE(half_space.is_inside(point));
  point.set_vector({-0.001, 0, 0});
  EXPECT_TRUE(half_space.is_inside(point));

  // finite size
  const double diameter = 1.;
  point.set_vector({-0.501, 0, 0});
  EXPECT_TRUE(half_space.is_inside(point, diameter));
  EXPECT_NEAR(half_space.nearest_distance(point), -0.501, NEAR_ZERO);
  point.set_vector({-0.499, 0, 0});
  EXPECT_FALSE(half_space.is_inside(point, diameter));
  EXPECT_NEAR(half_space.nearest_distance(point), -0.499, NEAR_ZERO);
//  point.set_vector({-0.501, 0, 0});
//  EXPECT_TRUE(half_space.is_inside(point, diameter));
//  EXPECT_NEAR(half_space.nearest_distance(point), -0.501, NEAR_ZERO);
}

TEST(HalfSpaceSine, amplitude1) {
  Position point;
  HalfSpaceSine half_space(
    MakeFormulaSineWave({{"width", "1"}}),
    { {"dimension", "0"},
      {"direction", "-1"},
      {"intersection", "0"},
      {"wave_dimension", "1"}});
  point.set_vector({0.001, 0, 0});
  EXPECT_FALSE(half_space.is_inside(point));
  point.set_vector({-0.001, 0, 0});
  EXPECT_TRUE(half_space.is_inside(point));
  point.set_vector({1.001, 0.25, 0});
  EXPECT_FALSE(half_space.is_inside(point));
  point.set_vector({0.999, 0.25, 0});
  EXPECT_TRUE(half_space.is_inside(point));
  point.set_vector({0.001, 0.5, 0});
  EXPECT_FALSE(half_space.is_inside(point));
  point.set_vector({-0.001, 0.5, 0});
  EXPECT_TRUE(half_space.is_inside(point));
  point.set_vector({-0.999, 0.75, 0});
  EXPECT_FALSE(half_space.is_inside(point));
  point.set_vector({-1.001, 0.75, 0});
  EXPECT_TRUE(half_space.is_inside(point));
  point.set_vector({0.001, 1, 0});
  EXPECT_FALSE(half_space.is_inside(point));
  point.set_vector({-0.001, 1, 0});
  EXPECT_TRUE(half_space.is_inside(point));

//  point.set_vector({1, 0, 0});
//  EXPECT_NEAR(half_space.nearest_distance(point), 0.4, 0.000001);
}

TEST(HalfSpaceSine, generalized) {
  Position point;
  HalfSpaceSine half_space(
    MakeFormulaSineWave({
      {"amplitude", "7.8"},
      {"width", "2400.3"}}),
    { {"dimension", "1"},
      {"direction", "1"},
      {"intersection", "1.3"},
      {"wave_dimension", "2"}});
  point.set_vector({0, 1.301, 0});
  EXPECT_TRUE(half_space.is_inside(point));
  point.set_vector({0, 1.299, 0});
  EXPECT_FALSE(half_space.is_inside(point));
  point.set_vector({0, 1.301+7.8, 2400.3/4.});
  EXPECT_TRUE(half_space.is_inside(point));
  point.set_vector({0, 1.299+7.8, 2400.3/4.});
  EXPECT_FALSE(half_space.is_inside(point));
  point.set_vector({0, 1.301-7.8, 3*2400.3/4.});
  EXPECT_TRUE(half_space.is_inside(point));
  point.set_vector({0, 1.299-7.8, 3*2400.3/4.});
  EXPECT_FALSE(half_space.is_inside(point));

  // finite size
  const double diameter = 1.;
  INFO("diameter " << diameter);
  point.set_vector({0, 1.301-7.8+0.5, 3*2400.3/4.});
  EXPECT_TRUE(half_space.is_inside(point, diameter));
  point.set_vector({0, 1.299-7.8+0.5, 3*2400.3/4.});
  EXPECT_FALSE(half_space.is_inside(point, diameter));

  // trough
  point.set_vector({0, 1.301+7.8+0.5, 2400.3/4.});
  EXPECT_TRUE(half_space.is_inside(point, diameter));
  point.set_vector({0, 1.299+7.8+0.5, 2400.3/4.});
  EXPECT_FALSE(half_space.is_inside(point, diameter));
}

TEST(HalfSpaceSine, integrate) {
  HalfSpaceSine half_space(
    MakeFormulaSineWave({
      {"amplitude", "1"},
      {"width", "1"}}),
    { {"dimension", "0"},
      {"intersection", "2"},
      {"direction", "-1"},
      {"wave_dimension", "1"}});
  auto random = MakeRandomMT19937();
  std::ofstream file("tmp/sineintxyz");
  for (double y = 0; y < 1; y += 0.1) {
    const double inte = half_space.integrate(
      Position({0.2, y, 0.}), random.get(), {
        {"wall_alpha0", "6"},
        {"wall_epsilon0", "-1"},
        {"wall_alpha1", "50"},
        {"wall_epsilon1", "0.01"},
        {"max_radius", "10"},
        {"num_shells", "100"},
        {"points_per_shell", "100"}});
    file << y << " " << inte << "\n";
  }
  //EXPECT_NEAR(inte, hamaker_half_plane(epsilon0, alpha0, distance), 0.04);
}

}  // namespace feasst
