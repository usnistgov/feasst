#include <vector>
#include "utils/test/utils.h"
#include "flat_histogram/include/collection_matrix.h"

namespace feasst {

// demonstrate that middle values in collection matrix affect ln_prob
TEST(CollectionMatrix, serialize) {
  TripleBandedCollectionMatrix colmat;
  colmat.resize(3);
  colmat.increment(0, 0, 0.2);
  colmat.increment(0, 1, 0.3);
  colmat.increment(0, 2, 0.3);
  colmat.increment(1, 0, 0.3);
  colmat.increment(1, 1, 0.3);
  colmat.increment(1, 2, 0.3);
  colmat.increment(2, 0, 0.3);
  colmat.increment(2, 1, 0.3);
  colmat.increment(2, 2, 0.3);
  DEBUG(feasst_str(colmat.matrix()));
  LnProbability lnpi;
  lnpi.resize(3);
  colmat.compute_ln_prob(&lnpi);
  DEBUG(feasst_str(lnpi.values()));
  colmat.increment(0, 1, 1000.3);
  colmat.increment(1, 1, 1000.3);
  colmat.increment(2, 1, 1000.3);
  DEBUG(feasst_str(colmat.matrix()));
  colmat.compute_ln_prob(&lnpi);
  DEBUG(feasst_str(lnpi.values()));
}

}  // namespace feasst
