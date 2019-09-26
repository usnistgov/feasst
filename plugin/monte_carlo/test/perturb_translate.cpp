#include "utils/test/utils.h"
#include "monte_carlo/include/perturb_translate.h"

namespace feasst {

TEST(PerturbTranslate, serialize) {
  PerturbTranslate add;
  PerturbTranslate add2 = test_serialize(add);
}

}  // namespace feasst
