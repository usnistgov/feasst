#include <sstream>
#include <gtest/gtest.h>
#include "core/include/log.h"

namespace feasst {

TEST(Log, serialize) {
  Log log;
  // log.set_steps_per(100);
  std::stringstream ss, ss2;
  log.serialize(ss);
  INFO(ss.str());
  // EXPECT_EQ("Log 497 0 0 -1 -1 0 1 0 0 1 ", ss.str());
  std::shared_ptr<Analyze> log2 = Log().deserialize(ss);
  log2->serialize(ss2);
  EXPECT_EQ(ss.str(), ss2.str());
}

}  // namespace feasst
