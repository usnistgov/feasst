
#ifndef FEASST_CORE_RANDOM_H_
#define FEASST_CORE_RANDOM_H_

#include <vector>
#include <random>
#include "core/include/position.h"
#include "core/include/debug.h"

namespace feasst {

/**
  Psuedo random number generator class.
  Initialize the random number seed by global function methods

  1. seed_random_by_date()
  2. seed_random(seed)
 */
class Random {
 public:
  Random() {
    generator_ = std::mt19937(rand());
    dis_double_ = std::uniform_real_distribution<double>(0.0, 1.0);
  }

  /// Return a random real number with a uniform probability distribution
  /// between 0 and 1.
  double uniform() {
    return dis_double_(generator_);
  }

  /// Return a random integer with a uniform probability distribution
  /// betwee min and max.
  int uniform(const int min, const int max);

  /// Return a random element within a vector.
  template<class T>
  T element(const std::vector<T> vector) {
    ASSERT(vector.size() > 0, "size error");
    return vector[uniform(0, vector.size() - 1)];
  }

  /// Return a random alpha numeric string of given length.
  std::string alpha_numeric(const int length = 5);

  /// Return a random position within a cube of side length with the origin
  /// at the center of the cube
  Position position_in_cube(const int dimension, const double length = 1.) {
    Position position;
    ASSERT(dimension == 3, "size error");
    position.set_to_origin_3D();
    for (int dim = 0; dim < dimension; ++dim) {
      position.set_coord(dim, (2.*uniform() - 1.)*length);
    }
    return position;
  }

  /// Return a random position within a cuboid of side lengths with the origin
  /// at the center
  Position position_in_cuboid(const Position side_length) {
    Position position;
    ASSERT(side_length.dimension() == 3, "size error");
    position.set_to_origin_3D();
    for (int dim = 0; dim < side_length.dimension(); ++dim) {
      position.set_coord(dim, (2.*uniform() - 1.)*side_length.coord(dim));
    }
    return position;
  }

 private:
  std::uniform_real_distribution<double> dis_double_;
  std::mt19937 generator_;
};

/// Initialize random number generator based on date and time.
void seed_random_by_date();

/// Initialize random number generator to seed value for reproducibility.
void seed_random(const int seed = 1346867550);

}  // namespace feasst

#endif  // FEASST_CORE_RANDOM_H_
