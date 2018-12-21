#include <gtest/gtest.h>
#include "core/test/particle_test.h"

namespace feasst {

TEST(Particle, getset) {
  Position position;
  Site site;
  std::vector<double> x = {3.5, 796.4, -45.4};
  position.set_vector(x);
  site.set_position(position);
}

TEST(Particle, check_size) {
  Particle particle = default_particle();
  Site site;
  Position pos;
  pos.set_vector({0, 0});
  site.set_position(pos);
  try {
    particle.add(site);
    particle.check_size();
    CATCH_PHRASE("size error");
  }
}

}  // namespace feasst
