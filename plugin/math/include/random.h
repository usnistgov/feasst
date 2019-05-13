
#ifndef FEASST_MATH_RANDOM_H_
#define FEASST_MATH_RANDOM_H_

#include <vector>
#include <random>
#include "math/include/position.h"
#include "utils/include/debug.h"
#include "math/include/constants.h"

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
    reseed_();
  }

  // copy constructor should re-seed
  Random(const Random& random) {
    reseed_();
  }

  /// Return a random real number with a uniform probability distribution
  /// between 0 and 1.
  double uniform() {
    return dis_double_(generator_);
  }

  /// Randomly return true or false
  bool coin_flip() {
    if (uniform() > 0.5) {
      return true;
    }
    return false;
  }

  /// Return a random integer with a uniform probability distribution
  /// betwee min and max.
  int uniform(const int min, const int max);

  /// Return random real number with a uniform probability distribution
  /// between min and max.
  double uniform_real(const double min, const double max);

  /// Return a random element within a vector.
  template<class T>
  T element(const std::vector<T> vector,
    /// optionally return the index associated with the element.
    int * index = NULL) {
    ASSERT(vector.size() > 0, "size error");
    const int vec_index = uniform(0, vector.size() - 1);
    if (index != NULL) {
      *index = vec_index;
    }
    return vector[vec_index];
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

  /// Random point on the surface of a unit sphere.
  void unit_sphere_surface(Position * position) {
    if (position->dimension() == 3) {
      // thanks to http://mathworld.wolfram.com/SpherePointPicking.html
      const double theta = 2*PI*uniform();
      const double phi = acos(2*uniform() - 1);
      PositionSpherical sph;
      sph.set_vector({1., theta, phi});
      *position = sph.cartesian();
    } else if (position->dimension() == 2) {
      const double theta = 2*PI*uniform();
      position->set_coord(0, cos(theta));
      position->set_coord(1, sin(theta));
    } else {
      ERROR("unrecognized dimension(" << position->dimension() << ")");
    }
  }

  /// Given a cumulative probability distribution, return a random integer
  /// index from a uniform probability distribution.
  /// The cumulative distribution must be monotonically nondecreasing.
  /// In addition, it must end with the value of unity.
  int index_from_cumulative_probability(const std::vector<double>& cumulative) {
    ASSERT(std::abs(cumulative.back() - 1) < NEAR_ZERO,
      "the cumulative distribution must end with a value of unity");
    const double random_uniform = uniform();
    double prev = -1;
    for (int index = 0; index < static_cast<int>(cumulative.size()); ++index) {
      const double prob = cumulative[index];
      ASSERT(prev <= prob, "cumulative probability must be monotonically " <<
        "nondecreasing");
      prev = prob;
      if (random_uniform < cumulative[index]) {
        return index;
      }
    }
    ERROR("should never reach this point");
    return -1;
  }

  void serialize(std::ostream& ostr) const {
    ostr << MAX_PRECISION;
    ostr << "1 "; // version
    ostr << generator_;
  }

  Random(std::istream& istr) {
    int version;
    istr >> version;
    istr >> generator_;
  }

 private:
  std::uniform_real_distribution<double> dis_double_;
  std::mt19937 generator_;
//  unsigned int seed_;

  // use by constructor and copy constructor
  void reseed_() {
    const int seed = rand();
    TRACE("seed " << seed << " address " << this);
    generator_ = std::mt19937(seed);
    dis_double_ = std::uniform_real_distribution<double>(0.0, 1.0);
  }
};

/// Initialize random number generator based on date and time.
void seed_random_by_date();

/// Initialize random number generator to seed value for reproducibility.
void seed_random(const int seed = 1346867550);

}  // namespace feasst

#endif  // FEASST_MATH_RANDOM_H_
