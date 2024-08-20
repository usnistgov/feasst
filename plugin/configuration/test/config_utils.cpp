#include "utils/include/io.h"
#include "utils/include/debug.h"
#include "configuration/test/config_utils.h"
#include "configuration/include/file_xyz.h"
#include "configuration/include/domain.h"

namespace feasst {

Configuration lj_sample4() {
  return *MakeConfiguration({{"particle_type0", "../particle/lj.fstprt"},
    {"xyz_file", "../plugin/configuration/test/data/lj_sample_config_periodic4.xyz"}});
}

Configuration spce_sample1() {
  return *MakeConfiguration({{"particle_type0", "../particle/spce.fstprt"},
    {"xyz_file", "../plugin/configuration/test/data/spce_sample_config_periodic1.xyz"}});
}

Configuration two_particle_configuration(argtype args) {
  const double cubic_side_length = dble("cubic_side_length", &args, 6.);
  Configuration config({{"cubic_side_length", feasst::str(cubic_side_length)},
                        {"particle_type0", "../particle/atom.fstprt"}});
  config.add_particle_of_type(0);
  config.add_particle_of_type(0);
  config.update_positions({{0, 0, 0}, {1.25, 0, 0}});
  feasst_check_all_used(args);
  return config;
}

}  // namespace feasst
