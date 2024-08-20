#include <cmath>  // abs
#include "utils/include/debug.h"
#include "utils/include/io.h"
#include "math/include/constants.h"
#include "configuration/test/config_utils.h"
#include "configuration/include/domain.h"
#include "system/test/sys_utils.h"
#include "system/include/system.h"
#include "system/include/lennard_jones.h"
#include "system/include/long_range_corrections.h"
#include "system/include/visit_model_cell.h"

namespace feasst {

System two_particle_system(argtype args) {
  System sys;
  sys.add(std::make_shared<Configuration>(two_particle_configuration(args)));
  sys.add_to_unoptimized(MakePotential(MakeLennardJones()));
  feasst_check_all_used(args);
  return sys;
}

}  // namespace feasst
