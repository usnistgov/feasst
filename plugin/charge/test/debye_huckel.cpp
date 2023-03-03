#include "utils/test/utils.h"
#include "charge/include/debye_huckel.h"

namespace feasst {

TEST(DebyeHuckel, serialize) {
  auto model = MakeDebyeHuckel({{"kappa", "0.5"}, {"dielectric", "1.2"}});
  std::shared_ptr<Model> model2 = test_serialize<DebyeHuckel, Model>(*model,
    "DebyeHuckel 2094 -1 -1 -1 -1 3682 0 0.5 1.2 -1 ");
}

}  // namespace feasst
