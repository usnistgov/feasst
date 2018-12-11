#include <gtest/gtest.h>
#include "core/include/position.h"
#include "core/include/site.h"
#include "core/include/particle.h"
#include "core/include/debug.h"

TEST(Particle, getset) {
  feasst::Position position;
  feasst::Site site;
  std::vector<double> x = {3.5, 796.4, -45.4};
  position.set_vector(x);
  site.set_position(position);
}

TEST(Particle, check_size) {
  try {
    feasst::Particle particle;
    particle.default_particle();
    feasst::Site site;
    feasst::Position pos;
    pos.set_vector({0, 0});
    site.set_position(pos);
    particle.add(site);
    particle.check_size();
    CATCH_PHRASE("size error");
  }
}
