#include <vector>
#include <memory>
#include "utils/test/utils.h"
#include "flat_histogram/include/collection_matrix.h"

namespace feasst {

// demonstrate that middle values in collection matrix affect ln_prob
TEST(CollectionMatrix, serialize) {
  CollectionMatrix colmat;
  colmat.resize(3);
  colmat.increment(0, 0, 0.2);
//  colmat.increment(0, 1, 0.3);
  colmat.increment(0, 1, 0.3);
  colmat.increment(1, 0, 0.3);
//  colmat.increment(1, 1, 0.3);
  colmat.increment(1, 1, 0.3);
  colmat.increment(2, 0, 0.3);
//  colmat.increment(2, 1, 0.3);
  colmat.increment(2, 1, 0.3);
  CollectionMatrix colmat2 = test_serialize(colmat);
  EXPECT_DOUBLE_EQ(colmat2.matrix()[0][0].sum(), 0.2);
  //DEBUG(feasst_str(colmat2.matrix()));
  LnProbability lnpi;
  lnpi.resize(3);
  colmat2.compute_ln_prob(&lnpi);
  DEBUG(feasst_str(lnpi.values()));
//  colmat2.increment(0, 1, 1000.3);
//  colmat2.increment(1, 1, 1000.3);
//  colmat2.increment(2, 1, 1000.3);
  CollectionMatrix colmat3 = test_serialize(colmat2);
  //DEBUG(feasst_str(colmat3.matrix()));
  colmat3.compute_ln_prob(&lnpi);
  DEBUG(feasst_str(lnpi.values()));
//  EXPECT_EQ(1, colmat3.min_blocks_());
}

//TEST(CollectionMatrix, blocks) {
//  auto cm = MakeCollectionMatrix();
//  cm->resize(6);
//  EXPECT_DOUBLE_EQ(cm->matrix()[0][0], 0);
//  EXPECT_DOUBLE_EQ(cm->matrix()[0][1], 0);
//  EXPECT_DOUBLE_EQ(cm->matrix()[0][2], 0);
//  cm->increment(0, 1, 1);
//  EXPECT_DOUBLE_EQ(cm->matrix()[0][1], 1.);
//  EXPECT_EQ(cm->blocks().size(), 5);
//  EXPECT_EQ(cm->blocks()[0].size(), 1);
//  EXPECT_DOUBLE_EQ(cm->blocks()[0][0].matrix()[0][1], 1.);
//  cm->increment(0, 1, 1);
//  cm->increment(0, 1, 1);
//  cm->increment(0, 1, 1);
//  cm->increment(0, 1, 1);
//  EXPECT_DOUBLE_EQ(cm->blocks()[0][0].matrix()[0][1], 1.);
//  EXPECT_DOUBLE_EQ(cm->blocks()[0][1].matrix()[0][1], 1.);
//  EXPECT_DOUBLE_EQ(cm->blocks()[0][4].matrix()[0][1], 1.);
//  EXPECT_EQ(cm->blocks()[0].size(), 5);
//  EXPECT_DOUBLE_EQ(cm->blocks()[1][0].matrix()[0][1], 2.);
//  EXPECT_DOUBLE_EQ(cm->blocks()[1][1].matrix()[0][1], 2.);
//  EXPECT_EQ(cm->blocks()[1].size(), 2);
//  EXPECT_DOUBLE_EQ(cm->blocks()[2][0].matrix()[0][1], 4.);
//  EXPECT_EQ(cm->blocks()[2].size(), 1);
//  EXPECT_EQ(cm->blocks()[3].size(), 0);
//  for (int i = 0; i < 1e2; ++i) cm->increment(0, 1, 1);
//  CollectionMatrix colmat2 = test_serialize(*cm);
//  EXPECT_EQ(colmat2.blocks()[3].size(), 3);
//}

}  // namespace feasst
