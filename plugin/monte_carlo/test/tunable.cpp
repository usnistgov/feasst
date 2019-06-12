#include "utils/test/utils.h"
#include "monte_carlo/include/tunable.h"

namespace feasst {

TEST(Tunable, serialize) {
  Tunable tunable;
  
  std::stringstream ss, ss2;
  tunable.serialize(ss);
  INFO(ss.str());
  Tunable tunable2(ss);
  tunable2.serialize(ss2);
  INFO(ss2.str());

  test_serialize(tunable);
}

}  // namespace feasst
