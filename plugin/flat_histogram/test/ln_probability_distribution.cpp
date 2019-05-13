#include "utils/test/utils.h"
#include "flat_histogram/include/ln_probability_distribution.h"

namespace feasst {

TEST(LnProbabilityDistribution, serialize) {
  LnProbabilityDistribution ln_prob;
  ln_prob.resize(5);
  ln_prob.add(0, 1.);
  LnProbabilityDistribution ln_prob2 = test_serialize(ln_prob, "885 5 1 0 0 0 0 ");
}

}  // namespace feasst
