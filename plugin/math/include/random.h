
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

  /// Return a random rotation matrix.
  RotationMatrix rotation(
    /// dimensionality of space
    const int dimension,
    /// maximum angle of rotation in degrees.
    const double max_angle = 180);

  /// Same as above, but optimized to reuse existing axis and matrix.
  void rotation(
    /// dimensionality of space
    const int dimension,
    Position * axis,
    RotationMatrix * rot_mat,
    /// maximum angle of rotation in degrees.
    const double max_angle = 180);

  /// Return the normal distribution with 0 mean and unit sigma (standard).
  /// Use Box-Muller transformation.
  /// See: https://mathworld.wolfram.com/Box-MullerTransformation.html
  virtual double standard_normal();

  /// Same as above, but with a given mean and standard deviation.
  double normal(const double mean, const double stdev) {
    return mean + stdev*standard_normal(); }

  /**
    Return a randomly selected bond length with harmonic potential of the form:
    betaU ~ spring_constant*(length - equilibrium_length)**2
    The typical 1/2 factor is included in spring_constant.
    prob(length) ~ length**2 exp(-betaU) dlength
    as described in Frenkel and Smit, Alg 43, page 578
    and Allen and Tildesley, Section G.3.
    The maximal length is 3 sigma beyond the mean.
    Only currently implemented for 3 dimensions. */
  double harmonic_bond_length(const double equilibrium_length,
    const double spring_constant,
    const int dimension);

  /**
    Same as above, but generalized for arbitrary exponential powers.
    betaU ~ spring_constant*(length - equilibrium_length)**exponent
    In this implementation, the maximum bond length is twice the equilibrium.
    If exponent == 2, use harmonic_bond_length. */
  double bond_length(const double equilibrium_length,
    const double spring_constant,
    const int exponent,
    const int dimension);

  /**
    Return bond angle selected from probability distribution associated with
    bending energy, \f$ \beta U=spring_constant*(t-equil_ang)^{exponent} \f$
    from Frenkel and Smit, page 343, below Equation 13.3.6.
    The equilibrium_angle has units of radians.
  */
  double bond_angle(const double equilibrium_angle,
    const double spring_constant,
    const int exponent,
    const int dimension,
    /// Optionally, disallow angles < minimum_angle. The max angle is PI.
    const double minimum_angle = 0.);

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
  std::map<std::string, std::shared_ptr<Random> >& deserialize_map();
  std::shared_ptr<Random> deserialize(std::istream& istr);
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
