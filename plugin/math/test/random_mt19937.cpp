#include <cmath>
#include <vector>
#include "utils/test/utils.h"
#include "math/include/random_mt19937.h"
#include "math/include/histogram.h"

namespace feasst {

TEST(Random, uniform) {
  RandomMT19937 random;
  random.seed_by_time();
  double sum = 0;
  const int num = 1e5;
  for (int i = 0; i < num; ++i) {
    sum += random.uniform();
  }
  EXPECT_NEAR(sum, double(num)/2., sqrt(double(num)));
}

TEST(Random, uniform_int) {
  RandomMT19937 random;
  random.seed_by_time();
  double sum = 0;
  const int num = 1e5;
  for (int i = 0; i < num; ++i) {
    sum += random.uniform(0, 10);
  }
  EXPECT_NEAR(sum, 5*double(num), 15*sqrt(double(num)));
}

TEST(Random, alpha_numeric) {
  RandomMT19937 random;
  random.seed_by_time();
  const int size = 10;
  std::string unique = random.alpha_numeric(size);
  INFO("unique alpha numeric: " << unique);
  EXPECT_EQ(unique.size(), size);
  EXPECT_NE(unique, random.alpha_numeric(size));
}

TEST(Random, element) {
  RandomMT19937 random;
  random.seed_by_time();
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
  random.seed_by_time();
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

TEST(Random, spherical_shell) {
  const double upper = 3, lower = 1.;
  RandomMT19937 random;
  random.seed_by_time();
  Position position;

  // 3D
  position.set_vector({0., 0., 0.});
  for (int point = 0; point < 1e3; ++point) {
    random.position_in_spherical_shell(lower, upper, &position);
    EXPECT_LE(lower, position.distance());
    EXPECT_GE(upper, position.distance());
    // visualize
    // std::cout << position.str() << std::endl;
  }

  // 2D
  position.set_vector({0., 0.});
  for (int point = 0; point < 1e3; ++point) {
    random.position_in_spherical_shell(lower, upper, &position);
    EXPECT_LE(lower, position.distance());
    EXPECT_GE(upper, position.distance());
    // visualize
    // std::cout << position.str() << std::endl;
  }
}

TEST(Random, index_from_cumulative_probability) {
  RandomMT19937 random;
  random.seed_by_time();
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
  random.seed_by_time();
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
  TRY(
    random3.set_cache_to_unload(random2);
    random3.uniform();
    random3.uniform();
    CATCH_PHRASE("can not unload if nothing stored");
  );
}

TEST(Random, standard_normal) {
  RandomMT19937 random;
  random.seed_by_time();
  double sum = 0., sumsq = 0.;
  const int num = 1e5;
  const double num_dble = static_cast<double>(num);
  for (int i = 0; i < num; ++i) {
    const double ran = random.standard_normal();
    sum += ran;
    sumsq += ran*ran;
  }
  const double av = sum/num_dble;
  const double stdev = std::sqrt((sumsq/num_dble-av*av));
  EXPECT_NEAR(sum, 0., 4*std::sqrt(num_dble));
  EXPECT_NEAR(stdev, 1., 0.01);
}

TEST(Random, bond_length) {
  RandomMT19937 random;
  random.seed_by_time();
  double spring_constant = 100.;
  const double equilibrum_length = 3.;
  double sum1 = 0., sumsq1 = 0., sum2 = 0., sumsq2 = 0.;
  const int num = 5e3;
  const double num_dble = static_cast<double>(num);
//  Histogram hist1, hist2;
//  hist1.set_width_center(0.01, equilibrum_length);
//  hist2.set_width_center(0.01, equilibrum_length);
  //({{"width", "0.1"}, {"max", "6"}});;
  for (int i = 0; i < num; ++i) {
    const double ran1 = random.bond_length(equilibrum_length, spring_constant, 2, 3);
    const double ran2 = random.harmonic_bond_length(equilibrum_length, spring_constant, 3);
//    hist1.add(ran1);
//    hist2.add(ran2);
    sum1 += ran1; sumsq1 += ran1*ran1; sum2 += ran2; sumsq2 += ran2*ran2;
  }
  const double av1 = sum1/num_dble;
  const double av2 = sum2/num_dble;
  const double stdev1 = std::sqrt((sumsq1/num_dble-av1*av1));
  const double stdev2 = std::sqrt((sumsq2/num_dble-av2*av2));
  const double stdev_exp = 1./std::sqrt(2*spring_constant);
  EXPECT_NEAR(av1, equilibrum_length, 8e-2);
  EXPECT_NEAR(av2, equilibrum_length, 8e-2);
  EXPECT_NEAR(stdev1, stdev_exp, 5e-2);
  EXPECT_NEAR(stdev2, stdev_exp, 5e-2);
//  for (int bin = 0; bin < hist.size(); ++bin) {
//    std::cout << hist.center_of_bin(bin) << " " << hist.histogram()[bin] << std::endl;
//  }
}

TEST(Random, bond_angle) {
  RandomMT19937 random;
  random.seed_by_time();
  double sum1 = 0., sumsq1 = 0., sum2 = 0., sumsq2 = 0.;
  const int num = 2e3;
  const double num_dble = static_cast<double>(num);
  Position point(3);
  for (int i = 0; i < num; ++i) {
    const double ran1 = random.bond_angle(0, 0, 2, 3);
    random.position_in_spherical_shell(0, 1., &point);
    const double ran2 = point.spherical().coord(2);
    sum1 += ran1; sumsq1 += ran1*ran1; sum2 += ran2; sumsq2 += ran2*ran2;
  }
  const double av1 = sum1/num_dble;
  const double av2 = sum2/num_dble;
  const double stdev1 = std::sqrt((sumsq1/num_dble-av1*av1));
  const double stdev2 = std::sqrt((sumsq2/num_dble-av2*av2));
  const double stdev_exp = 0.675;
  EXPECT_NEAR(av1, PI/2, 8e-2);
  EXPECT_NEAR(av2, PI/2, 8e-2);
  EXPECT_NEAR(stdev1, stdev_exp, 5e-2);
  EXPECT_NEAR(stdev2, stdev_exp, 5e-2);
}

}  // namespace feasst
