#include <sstream>
#include <gtest/gtest.h>
#include "core/include/criteria_writer.h"

namespace feasst {

TEST(CriteriaWriter, serialize) {
  CriteriaWriter analyze;
  std::stringstream ss;
  analyze.serialize(ss);
  EXPECT_EQ("CriteriaWriter 1 ", ss.str());
  std::shared_ptr<Analyze> analyze2 = CriteriaWriter().deserialize(ss);
  ss.str("");
  analyze2->serialize(ss);
  EXPECT_EQ("CriteriaWriter 1 ", ss.str());
}

}  // namespace feasst
