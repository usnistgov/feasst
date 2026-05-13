#include "utils/test/utils.h"
#include "configuration/include/configuration.h"
#include "aniso/include/recursive_table.h"

namespace feasst {

TEST(RecursiveTable, serialize) {
  //auto config = MakeConfiguration({{"particle_type", "atom:../particle/atom_new.txt"}});
  //const std::string table_file = "A_A:../plugin/aniso/test/data/dat_sqw_3rel_2z.txt";
  //auto vis = std::make_shared<RecursiveTable>(argtype({{"input_file", table_file}}));
  //vis->precompute(config.get());
  auto vis = std::make_shared<RecursiveTable>();
  auto vis2 = test_serialize(*vis);
}

}  // namespace feasst
