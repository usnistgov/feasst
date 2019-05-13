#include "utils/test/utils.h"
#include "system/include/system.h"
#include "system/test/system_test.h"

namespace feasst {

TEST(System, serialize) {
  System system = default_system();
  System system2 = test_serialize(system);
  EXPECT_EQ(system2.configuration().particle(1).site(0).position().coord(0), 1.25);
}

}  // namespace feasst
