#include <sstream>
#include <gtest/gtest.h>
#include "flat_histogram/include/ln_probability_distribution.h"
#include "debug.h"

namespace feasst {

TEST(LnProbabilityDistribution, serialize) {
  LnProbabilityDistribution ln_prob;
  ln_prob.resize(5);
  ln_prob.add(0, 1.);
  std::stringstream ss, ss2;
  ln_prob.serialize(ss);
  EXPECT_EQ("885 5 1 0 0 0 0 ", ss.str());
  LnProbabilityDistribution ln_prob2(ss);
  ln_prob2.serialize(ss2);
  EXPECT_EQ(ss.str(), ss2.str());
}

}  // namespace feasst
