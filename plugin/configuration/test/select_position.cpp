#include "utils/test/utils.h"
#include "configuration/include/select_position.h"
#include "configuration/include/configuration.h"

namespace feasst {

TEST(SelectPosition, serialize) {
  Configuration config;
  config.add_particle_type("../forcefield/data.spce");
  config.add_particle_of_type(0);
  SelectPosition sel(config.group_select(0), config.particles());
  SelectPosition sel2 = test_serialize(sel);
}

}  // namespace feasst
