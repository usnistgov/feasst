
#ifndef FEASST_SYSTEM_UTILS_H_
#define FEASST_SYSTEM_UTILS_H_

#include "utils/include/arguments.h"
#include "system/include/system.h"

namespace feasst {

/**
  Return a system with two Lennard-Jones particles.

  args:
  - cubic_box_length: default 6
 */
System two_particle_system(const argtype& args = argtype());

}  // namespace feasst

#endif  // FEASST_SYSTEM_UTILS_H_
