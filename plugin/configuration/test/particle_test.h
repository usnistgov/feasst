
#ifndef FEASST_CORE_PARTICLE_TEST_H_
#define FEASST_CORE_PARTICLE_TEST_H_

#include "configuration/include/particle.h"
#include "configuration/test/site_test.h"

namespace feasst {

inline Particle default_particle() {
  Particle part;
  Site site = default_site();
  part.add(site);
  Position pos = default_position();
  part.set_position(pos);
  part.set_type(0);
  return part;
}

}  // namespace feasst

#endif  // FEASST_CORE_PARTICLE_TEST_H_
