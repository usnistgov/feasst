#include <vector>
#include <memory>
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
  EXPECT_DOUBLE_EQ(colmat.matrix()[0][0], 0.2);
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

TEST(TripleBandedCollectionMatrix, blocks) {
  auto cm = MakeTripleBandedCollectionMatrix({{"min_block_size", "1"}});
  cm->resize(6);
  EXPECT_DOUBLE_EQ(cm->matrix()[0][0], 0);
  EXPECT_DOUBLE_EQ(cm->matrix()[0][1], 0);
  EXPECT_DOUBLE_EQ(cm->matrix()[0][2], 0);
  cm->increment(0, 1, 1);
  EXPECT_DOUBLE_EQ(cm->matrix()[0][1], 1.);
  EXPECT_EQ(cm->blocks().size(), 20);
  EXPECT_EQ(cm->blocks()[0].size(), 1);
  EXPECT_DOUBLE_EQ(cm->blocks()[0][0].matrix()[0][1], 1.);
  cm->increment(0, 1, 1);
  cm->increment(0, 1, 1);
  cm->increment(0, 1, 1);
  cm->increment(0, 1, 1);
  EXPECT_DOUBLE_EQ(cm->blocks()[0][0].matrix()[0][1], 1.);
  EXPECT_DOUBLE_EQ(cm->blocks()[0][1].matrix()[0][1], 1.);
  EXPECT_DOUBLE_EQ(cm->blocks()[0][4].matrix()[0][1], 1.);
  EXPECT_EQ(cm->blocks()[0].size(), 5);
  EXPECT_DOUBLE_EQ(cm->blocks()[1][0].matrix()[0][1], 2.);
  EXPECT_DOUBLE_EQ(cm->blocks()[1][1].matrix()[0][1], 2.);
  EXPECT_EQ(cm->blocks()[1].size(), 2);
  EXPECT_DOUBLE_EQ(cm->blocks()[2][0].matrix()[0][1], 4.);
  EXPECT_EQ(cm->blocks()[2].size(), 1);
  EXPECT_EQ(cm->blocks()[3].size(), 0);
}

}  // namespace feasst
