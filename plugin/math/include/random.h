
#ifndef FEASST_MATH_RANDOM_H_
#define FEASST_MATH_RANDOM_H_

#include <map>
#include <vector>
#include <deque>
#include <string>
#include <memory>
#include "utils/include/arguments.h"
#include "utils/include/cache.h"
#include "math/include/position.h"

namespace feasst {

class RotationMatrix;

/**
  Psuedo random number generator class.
  Note that all other Random distributions depend upon the uniform distribution,
  such that reproduction by storage is simplified.
 */
class Random {
 public:
  /**
    args:
    - seed : Provide an integer to seed the random number generator.
      If the string "time" is provided, then use the time and date to generate
      the seed.
      If no seed is provided, but random numbers are requested, then
      the time will be used to generate a seed.
      If the string "default" is provided, then use the default integer
      included in Random::seed().
   */
  explicit Random(argtype * args);

  /// Generate seed from time and date.
  void seed_by_time();

  /// Input seed by value for reproducibility.
  void seed(const int seed = 1346867550);

  /// Return a random real number with a uniform probability distribution
  /// between 0 and 1.
  double uniform();

  /// Return a random integer with a uniform probability distribution
  /// betwee min and max.
  int uniform(const int min, const int max);

  /// Randomly return true or false
  bool coin_flip();

  /// Return random real number with a uniform probability distribution
  /// between min and max.
  double uniform_real(const double min, const double max);

  /// Return a random index of a vector.
  template<class T>
  const int index(const std::vector<T>& vector) {
    // ASSERT(vector.size() > 0, "size error");
    const int vec_index = uniform(0, vector.size() - 1);
    return vec_index;
  }

  /// Return a pointer to a random element within a vector.
  template<class T>
  T * element(std::vector<T> * vector,
    /// optionally return the index associated with the element.
    int * return_index = NULL) {
    const int vec_index = index(*vector);
    if (return_index != NULL) {
      *return_index = vec_index;
    }
    return &(*vector)[vec_index];
  }

  /// Return a constant reference to a random element within a vector.
  template<class T>
  const T& const_element(const std::vector<T>& vector,
    /// optionally return the index associated with the element.
    int * return_index = NULL) {
    const int vec_index = index(vector);
    if (return_index != NULL) {
      *return_index = vec_index;
    }
    return vector[vec_index];
  }

  /// Return a random alpha numeric string of given length.
  std::string alpha_numeric(const int length = 5);

  /// Return a random position within a cube of side length with the origin
  /// at the center.
  Position position_in_cube(const int dimension, const double length = 1);

  /// Optimized version of the above, in that an existing position is modified.
  void position_in_cube(const int dimension,
                        const double length,
                        Position * position);

  /// Return a random position within a cuboid of side lengths with the origin
  /// at the center.
  Position position_in_cuboid(const Position& side_length);

  /// Optimized version of the above, in that an existing position is modified.
  void position_in_cuboid(const Position& side_length, Position * position);

  /// Random point on the surface of a unit sphere.
  void unit_sphere_surface(Position * position);

  // HWH rename position_in_shell (for 2D applications).
  /// Random point in a spherical shell.
  void position_in_spherical_shell(
    const double lower,
    const double upper,
    Position * position);

  /// Given a cumulative probability distribution, return a random integer
  /// index from a uniform probability distribution.
  /// The cumulative distribution must be monotonically nondecreasing.
  /// In addition, it must end with the value of unity.
  int index_from_cumulative_probability(const std::vector<double>& cumulative);

  /// Return a random quaternion using the method of:
  /// Franz J. Vesely, J. Comput. Phys., 47, 291-296 (1982).
  void quaternion(Position * quaternion);

  /// Return a random rotation matrix.
  RotationMatrix rotation(
    /// dimensionality of space
    const int dimension,
    /**
      If tunable == -1, generate a completely random rotation (default).
      Otherwise, In 3D, this is the relative weight of a random quaternion
      and a unit quaternion corresponding to no rotation.
      In 2D, this is the maximum rotation angle.
     */
    const double tunable = -1);

  /// Same as above, but optimized to reuse existing axis and matrix.
  void rotation(
    /// dimensionality of space
    const int dimension,
    Position * quaternion_or_axis,
    RotationMatrix * rot_mat,
    const double tunable = -1);

  /// Return the normal distribution with 0 mean and unit sigma (standard).
  /// Use Box-Muller transformation.
  /// See: https://mathworld.wolfram.com/Box-MullerTransformation.html
  virtual double standard_normal();

  /// Same as above, but with a given mean and standard deviation.
  double normal(const double mean, const double stdev) {
    return mean + stdev*standard_normal(); }

//  /**
//    Return bond angle selected from probability distribution associated with
//    bending energy as described in Frenkel and Smit, page 343,
//    below Equation 13.3.6.
//
//    \f$U(\theta)=\f$spring_constant\f$(\theta-\theta_0)^{exponent}\f$
//
//    \f$P(b)db \propto \exp[-\beta U(\theta)]d\cos\theta d\phi\f$
//
//    Note that the spring_constant has units of energy/degree^exponent.
//    The usual 1/2 factor is not included, but can be incorporated into
//    the spring constant manually by the user input to the particle file.
//  */
//  double bond_angle(const double theta0, /// units of radians
//    const double beta_spring_constant,  /// units of 1/degree^exponent
//    const int exponent,
//    const int dimension,
//    /// Optionally, disallow angles < minimum_angle. The max angle is PI.
//    const double minimum_angle = 0.);

  /// Return the cache.
  const Cache& cache() const { return cache_; }

  /// Set Cache to load.
  void set_cache_to_load(const bool load) { cache_.set_load(load); }

  /// Set Cache to unload.
  void set_cache_to_unload(const Random& random) {
    cache_.set_unload(random.cache()); }

  /// Serialize.
  std::string class_name() const { return class_name_; }
  virtual void serialize(std::ostream& ostr) const = 0;
  virtual std::shared_ptr<Random> create(std::istream& istr) const = 0;
  virtual std::shared_ptr<Random> create(argtype * args) const = 0;
  std::map<std::string, std::shared_ptr<Random> >& deserialize_map();
  std::shared_ptr<Random> deserialize(std::istream& istr);
  std::shared_ptr<Random> factory(const std::string name, argtype * args);
  virtual ~Random() {}

 protected:
  void parse_seed_(argtype * args);
  std::string class_name_ = "Random";
  void serialize_random_(std::ostream& ostr) const;
  explicit Random(std::istream& istr);

 private:
  Cache cache_;
  bool is_seeded_ = false;

  virtual void reseed_(const int seed) = 0;
  virtual double gen_uniform_() = 0;
  virtual int gen_uniform_(const int min, const int max);
};

}  // namespace feasst

#endif  // FEASST_MATH_RANDOM_H_
