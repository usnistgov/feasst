#include <iostream>
#include <fstream>
#include "utils/test/utils.h"
#include "utils/include/debug.h"
#include "math/include/random_mt19937.h"
#include "shape/include/formula_sine_wave.h"
#include "shape/include/slab_sine.h"

namespace feasst {

TEST(SlabSine, serialize) {
  auto slab = MakeSlabSine(MakeFormulaSineWave({{"amplitude", "2"}, {"width", "20"}}),
    { {"dimension", "0"}, {"wave_dimension", "1"}, {"average_bound0", "-5"},
      {"average_bound1", "5"}});
  std::shared_ptr<Shape> slab2 = test_serialize<SlabSine, Shape>(*slab);
  std::ofstream file("tmp/slabsinexyz");
  Position point;
  for (double x = -10; x <= 10; x += 0.5) {
  for (double y = -10; y <= 10; y += 0.5) {
  for (double z = -10; z <= 10; z += 0.5) {
    point.set_vector({x, y, z});
    if (slab2->is_inside(point)) {
      file << point.str() << std::endl;
    }
  }}}
}

TEST(SlabSine, integrate) {
  SlabSine slab(
    MakeFormulaSineWave({
      {"amplitude", "1"},
      {"width", "1"}}),
    { {"dimension", "0"},
      {"average_bound0", "-2"},
      {"average_bound1", "2"},
      {"wave_dimension", "1"}});
  auto random = MakeRandomMT19937();
  slab.integrate(
    Position({0.2, 0, 0.}), random.get(), {
      {"alpha0", "6"},
      {"epsilon0", "-1"},
      {"alpha1", "12"},
      {"epsilon1", "1"},
      {"max_radius", "10"},
      {"num_shells", "1000"},
      {"points_per_shell", "10000"}});
}

}  // namespace feasst
