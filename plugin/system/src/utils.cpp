#include "system/include/utils.h"
#include "configuration/include/utils.h"
#include "system/include/system.h"
#include "system/include/lennard_jones.h"

namespace feasst {

System two_particle_system(const argtype& args) {
  System sys;
  sys.add(two_particle_configuration(args));
  sys.add_to_unoptimized(Potential(MakeLennardJones()));
  return sys;
}

}  // namespace feasst
