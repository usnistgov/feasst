#include "utils/test/utils.h"
#include "system/include/synchronize_data.h"

namespace feasst {

TEST(SynchronizeData, serialize) {
  SynchronizeData data1;
  *data1.get_dble_1D() = {1, 2, 3};
  *data1.get_dble_2D() = {{1, 2, 3}, {4, 5.5, 6}};
  *data1.get_dble_5D() = {{{{{1, 2, 3}, {4, 5.5, 6}}}}};
  SynchronizeData data2 = test_serialize(data1);
  EXPECT_EQ(data2.dble_1D()[0], 1);
  EXPECT_EQ(data2.dble_1D()[1], 2);
  EXPECT_EQ(data2.dble_1D()[2], 3);
  EXPECT_EQ(data2.dble_2D()[0][0], 1);
  EXPECT_EQ(data2.dble_2D()[0][1], 2);
  EXPECT_EQ(data2.dble_2D()[0][2], 3);
  EXPECT_EQ(data2.dble_2D()[1][0], 4);
  EXPECT_EQ(data2.dble_2D()[1][1], 5.5);
  EXPECT_EQ(data2.dble_2D()[1][2], 6);
  EXPECT_EQ(data2.dble_5D()[0][0][0][0][0], 1);
  EXPECT_EQ(data2.dble_5D()[0][0][0][0][1], 2);
  EXPECT_EQ(data2.dble_5D()[0][0][0][0][2], 3);
  EXPECT_EQ(data2.dble_5D()[0][0][0][1][0], 4);
  EXPECT_EQ(data2.dble_5D()[0][0][0][1][1], 5.5);
  EXPECT_EQ(data2.dble_5D()[0][0][0][1][2], 6);
}

}  // namespace feasst
