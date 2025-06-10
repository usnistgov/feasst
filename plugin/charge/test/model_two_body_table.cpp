#include "utils/test/utils.h"
#include "math/include/table.h"
#include "configuration/include/configuration.h"
#include "configuration/include/domain.h"
#include "system/include/potential.h"
#include "system/include/model_two_body_table.h"
#include "system/include/lennard_jones.h"
#include "system/include/system.h"
#include "system/include/model_two_body_factory.h"
#include "charge/include/charge_screened.h"
#include "charge/include/coulomb.h"

namespace feasst {

TEST(ModelTwoBodyTable, rpm) {
  auto config = MakeConfiguration({{"cubic_side_length", "20"},
    {"particle_type", "rpm+:../plugin/charge/particle/rpm_plus.txt,rpm-:../plugin/charge/particle/rpm_minus.txt"},
    {"add_num_rpm+_particles", "1"},
    {"add_num_rpm-_particles", "1"}});
  config->update_positions({{0, 0, 0}, {0, 0, 2}});
  config->add_or_set_model_param("alpha", 0.2);
  auto model = MakeLennardJones();
  //auto model = MakeCoulomb();
  //auto model = MakeChargeScreened({{"table_size", "0"}});

  System no_table;
  no_table.add(config);
  no_table.add(MakePotential(model));
  no_table.precompute();

  System yes_table;
  yes_table.add(config);
  yes_table.add(MakePotential(model, {{"table_size", str(1e6)}}));
  yes_table.precompute();

  EXPECT_NEAR(no_table.energy(), yes_table.energy(), 1e-10);
}

}  // namespace feasst
