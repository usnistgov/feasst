#include <sstream>
#include <gtest/gtest.h>
#include "flat_histogram/include/bias_transition_matrix.h"

namespace feasst {

TEST(BiasTransitionMatrix, serialize) {
  auto bias = MakeBiasTransitionMatrix({{"min_sweeps", "1"}});
  bias->resize(Histogram({{"width", "1"}, {"max", "1"}}));
  bias->update(0, 0, 0, 0);
  std::stringstream ss, ss2;
  bias->serialize(ss);
  std::shared_ptr<Bias> bias2 = bias->deserialize(ss);
  bias2->serialize(ss2);
  EXPECT_EQ(ss.str(), ss2.str());
}

}  // namespace feasst
