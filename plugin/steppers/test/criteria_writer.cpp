#include <sstream>
#include <gtest/gtest.h>
#include "steppers/include/criteria_writer.h"

namespace feasst {

TEST(CriteriaWriter, serialize) {
  CriteriaWriter analyze;
  std::stringstream ss, ss2;
  analyze.serialize(ss);
  std::shared_ptr<Analyze> analyze2 = CriteriaWriter().deserialize(ss);
  analyze2->serialize(ss2);
  EXPECT_EQ(ss.str(), ss2.str());
}

}  // namespace feasst
