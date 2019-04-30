#include <gtest/gtest.h>
#include "core/include/system.h"
#include "core/test/system_test.h"

namespace feasst {

TEST(System, serialize) {
  System system = default_system();
  std::stringstream ss;
  system.serialize(ss);
  System system2(ss);
  EXPECT_EQ(system2.configuration().particle(1).site(0).position().coord(0), 1.25);
}

}  // namespace feasst
