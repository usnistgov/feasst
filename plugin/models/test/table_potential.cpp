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
  table->precompute(config.get());
  EXPECT_EQ(1, table->energy_table().size());
  EXPECT_EQ(1, table->energy_table()[0].size());
  auto table2 = test_serialize(*table);
  auto model = MakeTwoBodyTable();
  Position rel, pbc;
  rel.set_to_origin(3);
  pbc = rel;
  table2.compute(0, 0, 1, 0, config.get(), config->model_params(), model.get(), 1, &rel, &pbc);
  EXPECT_NEAR(-0.0615234375, table2.energy(), 1e-8);
}

}  // namespace feasst
