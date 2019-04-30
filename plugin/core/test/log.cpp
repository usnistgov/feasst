#include <sstream>
#include <gtest/gtest.h>
#include "core/include/log.h"

namespace feasst {

TEST(Log, serialize) {
  Log analyze;
  std::stringstream ss;
  analyze.serialize(ss);
  EXPECT_EQ("Log 1 ", ss.str());
  std::shared_ptr<Analyze> analyze2 = Log().deserialize(ss);
  ss.str("");
  analyze2->serialize(ss);
  EXPECT_EQ("Log 1 ", ss.str());
}

}  // namespace feasst
