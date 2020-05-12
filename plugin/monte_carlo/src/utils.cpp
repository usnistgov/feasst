#include "utils/include/utils_io.h"
#include "system/include/lennard_jones.h"
#include "system/include/long_range_corrections.h"
#include "monte_carlo/include/utils.h"
#include "monte_carlo/include/trial_add.h"
#include "monte_carlo/include/trial_remove.h"
#include "monte_carlo/include/metropolis.h"
#include "monte_carlo/include/trial_translate.h"
#include "monte_carlo/include/trial_rotate.h"

namespace feasst {

void add_trial_transfer(MonteCarlo * mc, const argtype& args) {
  mc->add(MakeTrialAdd(args));
  mc->add(MakeTrialRemove(args));
}

void lennard_jones(MonteCarlo * mc, const argtype& args) {
  Arguments args_(args);
  const double box_length = args_.key("cubic_box_length").dflt("8").dble();
  const std::string data = args_.key("particle").dflt("forcefield/data.lj").str();
  const bool translate = args_.key("translate").dflt("true").boolean();
  const std::string steps_per = args_.key("steps_per").dflt(str(1e4)).str();
  const bool lrc = args_.key("lrc").dflt("true").boolean();
//  const double cutoff = 2.;

  { System system;
    std::stringstream ss;
    ss << feasst::install_dir() << "/" << data;
    system.add(Configuration(
      MakeDomain({{"cubic_box_length", str(box_length)}}),
      {{"particle_type0", ss.str()},
      //{"init_cells", "1."}
    }));
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

    mc->set(system);
  }
  mc->set(MakeMetropolis({{"beta", "1.2"}, {"chemical_potential", "1."}}));
  if (translate) {
    mc->add(MakeTrialTranslate({{"weight", "1."}, {"tunable_param", "1."}}));
  }
}

}  // namespace feasst
