#include "utils/test/utils.h"
#include "math/include/table.h"
#include "configuration/include/configuration.h"
#include "configuration/include/domain.h"
#include "system/include/model_two_body_table.h"
#include "system/include/lennard_jones.h"
#include "system/include/system.h"
#include "system/include/model_two_body_factory.h"
#include "ewald/include/charge_screened.h"
#include "ewald/include/coulomb.h"

namespace feasst {

TEST(ModelTwoBodyTable, spce) {
  Configuration config(MakeDomain({{"cubic_box_length", "20"}}),
                       {{"particle_type0", "../plugin/ewald/forcefield/data.rpm_plus"},
                        {"particle_type1", "../plugin/ewald/forcefield/data.rpm_minus"}});
  config.add_particle_of_type(0);
  config.add_particle_of_type(1);
  config.update_positions({{0, 0, 0}, {0, 0, 2}});
  config.add_or_set_model_param("alpha", 0.2);
  //auto model = MakeLennardJones();
  auto model = MakeModelTwoBodyFactory({MakeLennardJones()});
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
