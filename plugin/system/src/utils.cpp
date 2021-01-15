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

System two_particle_system(const argtype& args) {
  System sys;
  sys.add(two_particle_configuration(args));
  sys.add_to_unoptimized(MakePotential(MakeLennardJones()));
  return sys;
}

System lennard_jones(const argtype& args) {
  Arguments args_(args);
  const double box_length = args_.key("cubic_box_length").dflt("8").dble();
  const std::string data = args_.key("particle").dflt("forcefield/data.lj").str();
  const bool lrc = args_.key("lrc").dflt("true").boolean();
  const double dual_cut = args_.key("dual_cut").dflt("-1").dble();
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
  return system;
}

}  // namespace feasst
