#include "utils/test/utils.h"
#include "chain/include/trials.h"

namespace feasst {

TEST(TrialGrow, serialize) {
  auto grow = MakeTrialGrow({
    {{"transfer", "true"},
     {"particle_type", "0"},
     {"site", "0"}}});
  Trial grow2 = test_serialize(*grow);
}

}  // namespace feasst
