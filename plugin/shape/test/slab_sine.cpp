#include <iostream>
#include <fstream>
#include "utils/test/utils.h"
#include "utils/include/debug.h"
#include "shape/include/formula_sine_wave.h"
#include "shape/include/slab_sine.h"

namespace feasst {

TEST(SlabSine, serialize) {
  SlabSine slab(MakeFormulaSineWave({{"amplitude", "2"}, {"width", "10"}, {"phase", "5"}}),
    { {"dimension", "0"}, {"wave_dimension", "1"}, {"average_bound0", "-5"},
      {"average_bound1", "5"}});
  SlabSine slab2 = test_serialize(slab);
  std::ofstream file("tmp/slabsinexyz");
  Position point;
  for (double x = -10; x <= 10; x += 0.5) {
  for (double y = -10; y <= 10; y += 0.5) {
  for (double z = -10; z <= 10; z += 0.5) {
    point.set_vector({x, y, z});
    if (slab2.is_inside(point)) {
      file << point.str() << std::endl;
    }
  }}}
}

}  // namespace feasst
