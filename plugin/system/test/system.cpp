#include "utils/test/utils.h"
#include "system/include/system.h"
#include "system/test/sys_utils.h"
#include "system/include/lennard_jones.h"
#include "system/include/long_range_corrections.h"

namespace feasst {

TEST(System, serialize) {
  System system = two_particle_system();
  System system2 = test_serialize(system);
  EXPECT_EQ(system2.configuration().particle(1).site(0).position().coord(0), 1.25);
}

TEST(System, cache) {
  System system = two_particle_system();
  System system2 = two_particle_system();
  DEBUG("enabling cache");
  system.load_cache(true);
  system.energy();
  DEBUG("cache loading? " << system.potential(0).cache().is_loading());
  DEBUG("cache loading? " << system.unoptimized().potential(0).cache().is_loading());
  system2.unload_cache(system);
  DEBUG(system.energy());
}

// compare with https://www.nist.gov/mml/csd/chemical-informatics-group/lennard-jones-fluid-reference-calculations-non-cuboid-cell
TEST(System, triclinic) {
  System system;
  system.add(MakeConfiguration({
    {"side_length0", "10"},
    {"side_length1", "9.84807753012208"},
    {"side_length2", "9.64974312607518"},
    {"xy", "1.7364817766693041"},
    {"xz", "2.5881904510252074"},
    {"yz", "0.42863479791864567"},
    {"particle_type0", "../forcefield/lj.fstprt"},
    {"xyz_file", "../plugin/configuration/test/data/lj_triclinic_sample_config_periodic3.xyz"}}));
  system.add(MakePotential(MakeLennardJones()));
  system.add(MakePotential(MakeLongRangeCorrections()));
  system.energy();
  EXPECT_DOUBLE_EQ(system.potential(0).stored_energy(), -505.78567945268367);
  EXPECT_DOUBLE_EQ(system.potential(1).stored_energy(), -29.37186430697248);
  INFO(system.stored_energy());
}

}  // namespace feasst
