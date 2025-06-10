#include "utils/test/utils.h"
#include "utils/include/debug.h"
#include "utils/include/debug.h"
#include "math/include/constants.h"
#include "math/include/utils_math.h"
#include "configuration/test/config_utils.h"
#include "configuration/include/file_xyz.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"

namespace feasst {

TEST(Configuration, tip4p) {
  auto config = std::make_unique<Configuration>(argtype({
    {"particle_type", "tip4p:../plugin/charge/particle/tip4p.txt"},
    {"add_num_tip4p_particles", "1"}
  }));
  EXPECT_EQ(1, config->num_particle_types());
  EXPECT_EQ("../plugin/charge/particle/tip4p.txt", config->type_to_file_name(0));
  EXPECT_EQ(3, config->num_site_types());
  EXPECT_EQ(0, config->site_type_to_particle_type(0));
  EXPECT_EQ(0, config->site_type_to_particle_type(1));
  EXPECT_EQ(0, config->site_type_to_particle_type(2));
  EXPECT_EQ("O", config->site_type_to_name(0));
  EXPECT_EQ("M", config->site_type_to_name(1));
  EXPECT_EQ("H", config->site_type_to_name(2));
  // INFO("config->str " << config->str());
  EXPECT_EQ(3, config->particle_type(0).num_bonds());
  EXPECT_EQ(2, config->unique_type(0).num_bonds());
  EXPECT_EQ(2, config->num_angle_types());
  EXPECT_EQ(3, config->particle_type(0).num_angles());
  EXPECT_EQ(2, config->unique_type(0).num_angles());
}

}  // namespace feasst
