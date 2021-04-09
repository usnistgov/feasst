#include <cmath>
#include "utils/test/utils.h"
#include "configuration/include/configuration.h"
#include "configuration/include/domain.h"
#include "system/include/lennard_jones.h"
#include "system/include/system.h"
#include "system/include/model_two_body_factory.h"
#include "system/include/model_two_body_table.h"

namespace feasst {

TEST(ModelTwoBodyTable, lj) {
  Configuration config(MakeDomain({{"cubic_box_length", "8"}}),
                       {{"particle_type0", "../forcefield/data.lj"},
                        {"particle_type1", "../forcefield/data.atom"}});
  config.add_particle_of_type(0);
  config.add_particle_of_type(1);
  config.update_positions({{0, 0, 0}, {0, 0, 2}});
  auto model = MakeLennardJones();
  //auto model = MakeModelTwoBodyFactory({MakeLennardJones()});
  System no_table;
  no_table.add(config);
  no_table.add(MakePotential(model));
  no_table.precompute();
  System yes_table;
  yes_table.add(config);
  yes_table.add(MakePotential(model, {{"table_size", str(1e3)}}));
  yes_table.precompute();
  EXPECT_NEAR(no_table.energy(), 4*(std::pow(2, -12) - std::pow(2, -6)), NEAR_ZERO);
  const double table_en = yes_table.energy();
  EXPECT_NE(table_en, 0);
  //EXPECT_NE(no_table.energy(), table_en);
  //EXPECT_NEAR(no_table.energy(), table_en, 1e-12);
}

}  // namespace feasst
