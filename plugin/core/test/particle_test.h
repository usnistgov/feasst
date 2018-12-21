
#ifndef FEASST_CORE_PARTICLE_TEST_H_
#define FEASST_CORE_PARTICLE_TEST_H_

#include "core/include/particle.h"
#include "core/test/site_test.h"

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
