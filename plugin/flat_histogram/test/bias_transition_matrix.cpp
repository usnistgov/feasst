#include "utils/test/utils.h"
#include "flat_histogram/include/bias_transition_matrix.h"

namespace feasst {

TEST(BiasTransitionMatrix, serialize) {
  auto bias = MakeBiasTransitionMatrix({{"min_sweeps", "1"}});
  bias->resize(Histogram({{"width", "1"}, {"max", "1"}}));
  bias->update(0, 0, 0, 0);
  std::shared_ptr<Bias> bias2 = test_serialize<BiasTransitionMatrix, Bias>(*bias);
}

}  // namespace feasst
