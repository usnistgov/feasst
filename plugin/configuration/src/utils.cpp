#include "utils/include/io.h"
#include "configuration/include/utils.h"
#include "configuration/include/file_xyz.h"
#include "configuration/include/domain.h"

namespace feasst {

Configuration lj_sample4() {
  auto config = MakeConfiguration({{"particle_type", "../forcefield/data.lj"}});
  FileXYZ().load(
    "../plugin/configuration/test/data/lj_sample_config_periodic4.xyz",
    config.get());
  return *config;
}

Configuration spce_sample1() {
  auto config = MakeConfiguration({{"particle_type", "../forcefield/data.spce"}});
  FileXYZ().load(
    "../plugin/configuration/test/data/spce_sample_config_periodic1.xyz",
    config.get());
  return *config;
}

Configuration two_particle_configuration(argtype args) {
  const double cubic_box_length = dble("cubic_box_length", &args, 6.);
  Configuration config(
    MakeDomain({{"cubic_box_length", feasst::str(cubic_box_length)}}),
    {{"particle_type", "../forcefield/data.atom"}});
  config.add_particle_of_type(0);
  config.add_particle_of_type(0);
  config.update_positions({{0, 0, 0}, {1.25, 0, 0}});
  check_all_used(args);
  return config;
}

}  // namespace feasst
