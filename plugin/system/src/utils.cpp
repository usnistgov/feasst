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
  sys.add_to_unoptimized(Potential(MakeLennardJones()));
  return sys;
}

System lennard_jones(const argtype& args) {
  Arguments args_(args);
  const double box_length = args_.key("cubic_box_length").dflt("8").dble();
  const std::string data = args_.key("particle").dflt("forcefield/data.lj").str();
  const bool lrc = args_.key("lrc").dflt("true").boolean();
  const double dual_cut = args_.key("dual_cut").dflt("-1").dble();
//  const double cutoff = 2.;

  System system;
  std::stringstream ss;
  ss << feasst::install_dir() << "/" << data;
  std::shared_ptr<Domain> domain;
  if (std::abs(dual_cut + 1) < NEAR_ZERO) {
    domain = MakeDomain({{"cubic_box_length", str(box_length)}});
  } else {
    domain = MakeDomain({{"cubic_box_length", str(box_length)},
                         {"init_cells", str(dual_cut)}});
  }
  system.add(Configuration(domain, {{"particle_type0", ss.str()}}));
  system.add(Potential(MakeLennardJones()));
//    { Potential potential(MakeLennardJones());
//      potential.set_model_params(system.configuration());
////      potential.set_model_param("cutoff", 0, cutoff);
////      EXPECT_NEAR(potential.model_params().mixed_cutoff()[0][0], cutoff, NEAR_ZERO);
//      system.add_to_unoptimized(potential); }
  if (lrc) {
    system.add(Potential(MakeLongRangeCorrections()));
  }
//    { Potential lrc(MakeLongRangeCorrections());
//      lrc.set_model_params(system.configuration());
////      lrc.set_model_param("cutoff", 0, cutoff);
////      EXPECT_NEAR(lrc.model_params().mixed_cutoff()[0][0], cutoff, NEAR_ZERO);
//      //system.add_to_unoptimized(lrc);
//    }
  if (std::abs(dual_cut + 1) > NEAR_ZERO) {
    Potential ref(MakeLennardJones(), MakeVisitModelCell());
    ref.set_model_params(system.configuration());
    ref.set_model_param("cutoff", 0, dual_cut);
    system.add_to_reference(ref);
  }
  return system;
}

}  // namespace feasst
