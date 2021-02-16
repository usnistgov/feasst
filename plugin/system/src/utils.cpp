#include <cmath>  // abs
#include "utils/include/debug.h"
#include "utils/include/io.h"
#include "math/include/constants.h"
#include "configuration/include/utils.h"
#include "configuration/include/domain.h"
#include "system/include/utils.h"
#include "system/include/system.h"
#include "system/include/lennard_jones.h"
#include "system/include/long_range_corrections.h"
#include "system/include/visit_model_cell.h"

namespace feasst {

System two_particle_system(argtype args) {
  System sys;
  sys.add(two_particle_configuration(args));
  sys.add_to_unoptimized(MakePotential(MakeLennardJones()));
  check_all_used(args);
  return sys;
}

System lennard_jones(argtype args) {
  const double box_length = dble("cubic_box_length", &args, 8);
  const std::string data = str("particle", &args, "forcefield/data.lj");
  const bool lrc = boolean("lrc", &args, true);
  const double dual_cut = dble("dual_cut", &args, -1.);
  std::shared_ptr<Domain> domain;
  System system;
  std::stringstream ss;
  ss << feasst::install_dir() << "/" << data;
  system.add(Configuration(MakeDomain({{"cubic_box_length", str(box_length)}}),
                           {{"particle_type0", ss.str()}}));
  system.add(MakePotential(MakeLennardJones()));
  if (lrc) system.add(MakePotential(MakeLongRangeCorrections()));
  if (std::abs(dual_cut + 1) > NEAR_ZERO) {
    auto ref = MakePotential(MakeLennardJones(),
                             MakeVisitModelCell({{"min_length", str(dual_cut)}}));
    ref->set_model_params(system.configuration());
    ref->set_model_param("cutoff", 0, dual_cut);
    system.add_to_reference(ref);
  }
  check_all_used(args);
  return system;
}

}  // namespace feasst
