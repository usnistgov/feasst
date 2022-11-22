#include <cmath>
#include "utils/test/utils.h"
#include "math/include/table.h"
#include "math/include/constants.h"

namespace feasst {

TEST(Table1D, interpolate) {
  auto table = MakeTable1D({{"num", "2"}, {"default_value", "0."}});
  EXPECT_NEAR(table->linear_interpolation(0.5), 0, NEAR_ZERO);
  table->set_data(0, 1.);
  auto table2 = std::make_shared<Table1D>(test_serialize(*table));
  EXPECT_NEAR(table2->linear_interpolation(0.5), 1./2., NEAR_ZERO);
  EXPECT_EQ(table2->bin_to_value(1), 1);
  EXPECT_EQ(table2->value_to_nearest_bin(0.45), 0);
  EXPECT_EQ(table2->value_to_nearest_bin(0.501), 1);
}

TEST(Table1D, forward_difference_interpolate) {
  auto table = MakeTable1D({{"num", str(11)}, {"default_value", "0."}});
  for (int bin = 0; bin < table->num(); ++bin) {
    const double x = table->bin_to_value(bin);
    table->set_data(bin, std::erfc(x));
  }
  Table1D table2 = test_serialize(*table);
  for (double x = 0; x <= 0.9; x += 0.001) {
    const double diff_lin = std::erfc(x) - table->linear_interpolation(x);
    const double diff_fd = std::erfc(x) - table->forward_difference_interpolation(x);

    EXPECT_LT(std::abs(diff_lin), 0.002);
    EXPECT_LT(std::abs(diff_fd ), 0.0005);
    //INFO(x << " " << std::erfc(x) << " " << diff_lin << " " << diff_fd);
  }
}

TEST(Table2D, interpolate) {
  auto table = MakeTable2D(
    {{"num0", "2"}, {"num1", "2"}, {"default_value", "0."}});
  EXPECT_NEAR(table->linear_interpolation(0.5, 0.5), 0, NEAR_ZERO);
  table->set_data(0, 0, 1.);
  auto table2 = std::make_shared<Table2D>(test_serialize(*table));
  EXPECT_NEAR(table2->linear_interpolation(0.5, 0.5), 1./4., NEAR_ZERO);
  EXPECT_EQ(table2->bin_to_value(0, 1), 1);
  EXPECT_EQ(table2->value_to_nearest_bin(0, 0.45), 0);
  EXPECT_EQ(table2->value_to_nearest_bin(0, 0.501), 1);
}

TEST(Table3D, interpolate) {
  auto table = MakeTable3D(
    {{"num0", "2"}, {"num1", "2"}, {"num2", "2"}, {"default_value", "0."}});
  EXPECT_NEAR(table->linear_interpolation(0.5, 0.5, 0.5), 0, NEAR_ZERO);
  table->set_data(0, 0, 0, 1.);
  auto table2 = std::make_shared<Table3D>(test_serialize(*table));
  EXPECT_NEAR(table2->linear_interpolation(0.5, 0.5, 0.5), 1./8., NEAR_ZERO);
  EXPECT_EQ(table2->bin_to_value(0, 1), 1);
  EXPECT_EQ(table2->value_to_nearest_bin(0, 0.45), 0);
  EXPECT_EQ(table2->value_to_nearest_bin(0, 0.501), 1);
  table2->write("tmp/tab3dint");
  auto table3 = MakeTable3D("tmp/tab3dint");
  EXPECT_NEAR(table3->linear_interpolation(0.5, 0.5, 0.5), 1./8., NEAR_ZERO);
  EXPECT_EQ(table3->bin_to_value(0, 1), 1);
  EXPECT_EQ(table3->value_to_nearest_bin(0, 0.45), 0);
  EXPECT_EQ(table3->value_to_nearest_bin(0, 0.501), 1);
}

TEST(Table4D, interpolate) {
  auto table = MakeTable4D({{"num0", "2"}, {"num1", "2"}, {"num2", "2"},
    {"num3", "2"}, {"default_value", "0."}});
  EXPECT_NEAR(table->linear_interpolation(0.5, 0.5, 0.5, 0.5), 0, NEAR_ZERO);
  table->set_data(0, 0, 0, 0, 1.);
  auto table2 = std::make_shared<Table4D>(test_serialize(*table));
  EXPECT_NEAR(table2->linear_interpolation(0.5, 0.5, 0.5, 0.5), 1./std::pow(2, 4), NEAR_ZERO);
  EXPECT_EQ(table2->bin_to_value(0, 1), 1);
  EXPECT_EQ(table2->value_to_nearest_bin(0, 0.45), 0);
  EXPECT_EQ(table2->value_to_nearest_bin(0, 0.501), 1);
}

TEST(Table5D, interpolate) {
  auto table = MakeTable5D({{"num0", "2"}, {"num1", "2"}, {"num2", "2"},
    {"num3", "2"}, {"num4", "2"}, {"default_value", "0."}});
  EXPECT_NEAR(table->linear_interpolation(0.5, 0.5, 0.5, 0.5, 0.5), 0, NEAR_ZERO);
  table->set_data(0, 0, 0, 0, 0, 1.);
  auto table2 = std::make_shared<Table5D>(test_serialize(*table));
  EXPECT_NEAR(table2->linear_interpolation(0.5, 0.5, 0.5, 0.5, 0.5), 1./std::pow(2, 5), NEAR_ZERO);
  EXPECT_EQ(table2->bin_to_value(0, 1), 1);
  EXPECT_EQ(table2->value_to_nearest_bin(0, 0.45), 0);
  EXPECT_EQ(table2->value_to_nearest_bin(0, 0.501), 1);
}

TEST(Table6D, interpolate) {
  auto table = MakeTable6D({{"num0", "2"}, {"num1", "2"}, {"num2", "2"},
    {"num3", "2"}, {"num4", "2"}, {"num5", "2"}, {"default_value", "0."}});
  EXPECT_NEAR(table->linear_interpolation(0.5, 0.5, 0.5, 0.5, 0.5, 0.5), 0, NEAR_ZERO);
  table->set_data(0, 0, 0, 0, 0, 0, 1.);
  auto table2 = std::make_shared<Table6D>(test_serialize(*table));
  EXPECT_NEAR(table2->linear_interpolation(0.5, 0.5, 0.5, 0.5, 0.5, 0.5), 1./std::pow(2, 6), NEAR_ZERO);
  EXPECT_EQ(table2->bin_to_value(0, 1), 1);
  EXPECT_EQ(table2->value_to_nearest_bin(0, 0.45), 0);
  EXPECT_EQ(table2->value_to_nearest_bin(0, 0.501), 1);
}

}  // namespace feasst
