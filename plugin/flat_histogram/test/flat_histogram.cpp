#include "utils/test/utils.h"
#include "flat_histogram/include/flat_histogram.h"
#include "flat_histogram/test/flat_histogram_test.h"

namespace feasst {

TEST(FlatHistogram, reweight) {
  std::shared_ptr<FlatHistogram> crit = crit_fh(0);
  LnProbability lnpirw = crit->reweight(1.5);
  EXPECT_NEAR(-7.75236, lnpirw.values()[0], 1e-5);
  test_serialize(*crit);
}

}  // namespace feasst
