#include "utils/test/utils.h"
#include "system/include/system.h"
#include "system/test/system_test.h"

namespace feasst {

TEST(System, serialize) {
  System system = default_system();
  System system2 = test_serialize(system);
  EXPECT_EQ(system2.configuration().particle(1).site(0).position().coord(0), 1.25);
}

TEST(System, cache) {
  System system = default_system();
  System system2 = default_system();
  DEBUG("enabling cache");
  system.load_cache(true);
  system.energy();
  DEBUG("cache loading? " << system.potential(0).cache().is_loading());
  DEBUG("cache loading? " << system.unoptimized().potentials()[0].cache().is_loading());
  system2.unload_cache(system);
  DEBUG(system.energy());
}

}  // namespace feasst
