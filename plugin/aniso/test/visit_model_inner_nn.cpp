#include "utils/test/utils.h"
#include "configuration/include/configuration.h"
#include "aniso/include/visit_model_inner_nn.h"

namespace feasst {

TEST(VisitModelInnerNN, serialize) {
  auto config = MakeConfiguration({{"particle_type0", "../particle/atom.txt"}});
  const std::string table_file = "../plugin/aniso/test/data/dat_sqw_3rel_2z.txt";
  auto vis = std::make_shared<VisitModelInnerNN>(argtype({{"table_file", table_file}}));
  vis->precompute(config.get());
  auto vis2 = test_serialize(*vis);
}

}  // namespace feasst
