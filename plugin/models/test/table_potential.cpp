#include "utils/test/utils.h"
#include "configuration/include/configuration.h"
#include "models/include/two_body_table.h"
#include "models/include/table_potential.h"

namespace feasst {

TEST(TablePotential, serialize) {
  auto config = MakeConfiguration({
    {"cubic_side_length", "8"},
    {"particle_type0", "../particle/lj.fstprt"},
    {"add_particles_of_type0", "2"}});
  config->update_positions({{0, 0, 0}, {2, 0, 0}});
  auto table = MakeTablePotential({{"table_file", "../plugin/models/test/data/lj_table.txt"}});
  table->precompute(config->model_params());
  TablePotential table2 = test_serialize(*table);
  EXPECT_EQ(1, table2.energy_table().size());
  EXPECT_EQ(1, table2.energy_table()[0].size());
  EXPECT_NEAR(-0.0615234375, table2.energy(4, 0, 0, config->model_params()), 1e-8);
}

}  // namespace feasst
