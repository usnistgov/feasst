#include "utils/test/utils.h"
#include "flat_histogram/include/transition_matrix.h"

namespace feasst {

TEST(TransitionMatrix, serialize) {
  auto bias = MakeTransitionMatrix({{"min_sweeps", "1"}});
//  bias->resize(Histogram({{"width", "1"}, {"max", "1"}}));
//  bias->update(0, 0, 0., false, true);
  auto bias2 = test_serialize_unique(*bias);
}

TEST(TransitionMatrix, args) {
  TRY(
    MakeTransitionMatrix();
    CATCH_PHRASE("key(min_sweeps) is required");
  );
}

}  // namespace feasst
