
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
System two_particle_system(argtype args = argtype());

/**
  Return a Lennard-Jones system

  args:
  - cubic_box_length: default 8
  - particle: default forcefield/data.lj
  - lrc: use long range corrections (default: true)
  - dual_cut: Add a reference potential using this short cutoff with a cell
    list.
    If -1, ignore. (default: -1).
 */
System lennard_jones(argtype args = argtype());

}  // namespace feasst

#endif  // FEASST_SYSTEM_UTILS_H_
