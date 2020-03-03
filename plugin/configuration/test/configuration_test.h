
#ifndef FEASST_CORE_CONFIGURATION_TEST_H_
#define FEASST_CORE_CONFIGURATION_TEST_H_

#include <vector>
#include <utility>
#include <string>
#include "configuration/include/configuration.h"
#include "configuration/include/file_xyz.h"

namespace feasst {

inline Configuration default_configuration() {
  auto config = MakeConfiguration(
    MakeDomain({{"cubic_box_length", "6"}}),
    {{"particle_type", "../forcefield/data.atom"}}
  );
  config->add_particle_of_type(0);
  config->add_particle_of_type(0);
  Position position = config->particle(1).site(0).position();
  position.set_coord(0, 1.25);
  Particle part1 = config->particle(1);
  part1.displace(position);
  Select select;
  select.add_sites(1, {0});
  config->replace_position(select, part1);
  return *config;
}

inline Configuration lj_sample() {
  auto config = MakeConfiguration({{"particle_type", "../forcefield/data.lj"}});
  FileXYZ().load("../plugin/system/test/data/lj_sample_config_periodic4.xyz",
    config.get());
  return *config;
}

inline Configuration spce_sample() {
  auto config = MakeConfiguration({{"particle_type", "../forcefield/data.spce"}});
  FileXYZ().load("../plugin/system/test/data/spce_sample_config_periodic1.xyz", config.get());
  return *config;
}

inline Configuration chain10_sample() {
  auto config = MakeConfiguration(
    MakeDomain({{"cubic_box_length", "10"}}),
    {{"particle_type", "../forcefield/data.chain10"}});
  config->add_particle_of_type(0);
  return *config;
}

}  // namespace feasst

#endif  // FEASST_CORE_CONFIGURATION_TEST_H_
