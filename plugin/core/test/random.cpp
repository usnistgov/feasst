#include <gtest/gtest.h>
#include "core/include/random.h"

namespace feasst {

TEST(Random, uniform) {
  Random random;
  double sum = 0;
  const int num = 1e5;
  for (int i = 0; i < num; ++i) {
    sum += random.uniform();
  }
  EXPECT_NEAR(sum, double(num)/2., sqrt(double(num)));
}

TEST(Random, uniform_int) {
  Random random;
  double sum = 0;
  const int num = 1e5;
  for (int i = 0; i < num; ++i) {
    sum += random.uniform(0, 10);
  }
  EXPECT_NEAR(sum, 5*double(num), 15*sqrt(double(num)));
}

TEST(Random, alpha_numeric) {
  Random random;
  const int size = 10;
  std::string unique = random.alpha_numeric(size);
  INFO("unique alpha numeric: " << unique);
  EXPECT_EQ(unique.size(), size);
  EXPECT_NE(unique, random.alpha_numeric(size));
}

TEST(Random, element) {
  Random random;
  std::vector<int> vec = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  double sum = 0;
  const int num = 1e5;
  for (int i = 0; i < num; ++i) {
    sum += random.element(vec);
  }
  EXPECT_NEAR(sum, 5*double(num), 15*sqrt(double(num)));
}

}  // namespace feasst
