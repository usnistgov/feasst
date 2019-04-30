#include <sstream>
#include <gtest/gtest.h>
#include "core/include/potential.h"

namespace feasst {

TEST(Potential, serialize) {
  Potential potential;
  std::stringstream ss, ss2;
  potential.serialize(ss);
  Potential potential2(ss);
  potential2.serialize(ss2);
  EXPECT_EQ(ss.str(), ss2.str());
}

}  // namespace feasst
