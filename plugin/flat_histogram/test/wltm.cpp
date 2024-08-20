#include "utils/test/utils.h"
#include "flat_histogram/include/wang_landau.h"
#include "flat_histogram/include/transition_matrix.h"
#include "flat_histogram/include/wltm.h"

namespace feasst {

TEST(WLTM, serialize) {
  auto bias = MakeWLTM({{"collect_flatness", "1"},
                        {"min_flatness", "2"},
                        {"min_sweeps", "3"}});
//  bias->resize(Histogram({{"width", "1"}, {"max", "1"}}));
//  bias->update(0, 0, 0., false, true);
  auto bias2 = test_serialize_unique(*bias);
}

}  // namespace feasst
