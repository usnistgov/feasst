#include "utils/test/utils.h"
#include "flat_histogram/include/wltm.h"

namespace feasst {

TEST(WLTM, serialize) {
  auto bias = MakeWLTM({{"collect_flatness", "1"},
                        {"min_flatness", "2"},
                        {"min_sweeps", "3"}});
  bias->resize(Histogram({{"width", "1"}, {"max", "1"}}));
  bias->update(0, 0, 0., false);
  std::shared_ptr<Bias> bias2 = test_serialize<WLTM, Bias>(*bias);
}

}  // namespace feasst
