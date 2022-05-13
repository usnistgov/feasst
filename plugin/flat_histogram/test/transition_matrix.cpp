#include "utils/test/utils.h"
#include "flat_histogram/include/transition_matrix.h"

namespace feasst {

TEST(TransitionMatrix, serialize) {
  auto bias = MakeTransitionMatrix({{"min_sweeps", "1"}});
//  bias->resize(Histogram({{"width", "1"}, {"max", "1"}}));
//  bias->update(0, 0, 0., false, true);
  std::shared_ptr<Bias> bias2 = test_serialize<TransitionMatrix, Bias>(*bias);
}

TEST(TransitionMatrix, args) {
  TRY(
    MakeTransitionMatrix();
    CATCH_PHRASE("key(min_sweeps) is required");
  );
}

}  // namespace feasst
