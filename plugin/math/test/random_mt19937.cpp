#include <vector>
#include <gtest/gtest.h>
#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"

namespace feasst {

TEST(Random, uniform) {
  RandomMT19937 random;
  double sum = 0;
  const int num = 1e5;
  for (int i = 0; i < num; ++i) {
    sum += random.uniform();
  }
  EXPECT_NEAR(sum, double(num)/2., sqrt(double(num)));
}

TEST(Random, uniform_int) {
  RandomMT19937 random;
  double sum = 0;
  const int num = 1e5;
  for (int i = 0; i < num; ++i) {
    sum += random.uniform(0, 10);
  }
  EXPECT_NEAR(sum, 5*double(num), 15*sqrt(double(num)));
}

TEST(Random, alpha_numeric) {
  RandomMT19937 random;
  const int size = 10;
  std::string unique = random.alpha_numeric(size);
  INFO("unique alpha numeric: " << unique);
  EXPECT_EQ(unique.size(), size);
  EXPECT_NE(unique, random.alpha_numeric(size));
}

TEST(Random, element) {
  RandomMT19937 random;
  std::vector<int> vec = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  double sum = 0;
  const int num = 1e5;
  for (int i = 0; i < num; ++i) {
    sum += random.const_element(vec);
  }
  EXPECT_NEAR(sum, 5*double(num), 15*sqrt(double(num)));
}

TEST(Random, unit_sphere) {
  RandomMT19937 random;
  Position position;
  position.set_vector({0., 0.});
  random.unit_sphere_surface(&position);
  EXPECT_NEAR(position.distance(), 1., NEAR_ZERO);
  position.set_vector({0., 0., 0.});
  random.unit_sphere_surface(&position);
  EXPECT_NEAR(position.distance(), 1., NEAR_ZERO);

//  // visualize
//  for (int point = 0; point < 1e3; ++point) {
//    random.unit_sphere_surface(&position);
//    std::cout << position.str() << std::endl;
//  }
}

TEST(Random, index_from_cumulative_probability) {
  RandomMT19937 random;
  std::vector<double> cpdf;
  const int ncpdf = 10, num = 100;
  for (int i = 0; i < ncpdf; ++i) cpdf.push_back((i+1)/double(ncpdf));
  std::vector<double> cpdfran(ncpdf);
  for (int i = 0; i < num; ++i) {
    const int j = random.index_from_cumulative_probability(cpdf);
    ++cpdfran[j];
  }
  for (int i = 0; i < ncpdf; ++i) EXPECT_NEAR(cpdfran[i]/double(num), ncpdf/double(num), 0.2);
}

TEST(Random, serialize) {
//  RandomMT19937 random;
//  std::shared_ptr<Random> random2 = test_serialize(random);
  RandomMT19937 random;
  random.set_cache_to_load(true);
  RandomMT19937 random2 = test_serialize(random);
  // EXPECT_NEAR(0.99971855963527156, random.uniform(), 1e-6);
  const double next = random.uniform();
  // EXPECT_NEAR(0.80124978724913187, next, 1e-6);
  EXPECT_EQ(next, random2.uniform());
  RandomMT19937 random3;
  random3.set_cache_to_unload(random2);
  EXPECT_EQ(next, random3.uniform());
  random3.set_cache_to_load(false);
  EXPECT_NE(random.uniform(), random3.uniform());
  try {
    random3.set_cache_to_unload(random2);
    random3.uniform();
    random3.uniform();
    CATCH_PHRASE("can not unload if nothing stored");
  }
}

}  // namespace feasst
