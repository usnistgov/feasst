#include <cmath>  // acos
#include <iostream>
#include "utils/include/arguments.h"
#include "utils/include/io.h"
#include "utils/include/serialize_extra.h"
#include "utils/include/cache.h"
#include "math/include/position.h"
#include "math/include/utils_math.h"
#include "math/include/constants.h"
#include "math/include/matrix.h"
#include "math/include/random.h"

namespace feasst {

Random::Random(argtype * args) {
  cache_ = std::make_shared<Cache>();
}

// parsing seed in constructor leads to pure virtual function reseed_
void Random::parse_seed_(argtype * args) {
  if (used("seed", *args)) {
    const std::string seed_str = str("seed", args);
    if (seed_str == "time") {
      WARN("Random argument \"time\" is deprecated and will be removed.");
      seed_by_time();
    } else if (seed_str == "default") {
      WARN("Random argument \"default\" is deprecated and will be removed.");
      seed();
    } else {
      seed(str_to_int(seed_str));
    }
  }
}

void Random::seed_by_time() {
  const int t = time(NULL);
  srand(t);
  std::cout << "# initializing random number generator with seed: "
            << t << std::endl;
  reseed_(t);
  is_seeded_ = true;
}

void Random::seed(const int seed) {
  srand(seed);
  std::cout << "# initializing random number generator with seed: "
            << seed << std::endl;
  reseed_(seed);
  is_seeded_ = true;
}

double Random::uniform() {
  if (!is_seeded_) {
    seed_by_time();
  }
  double ran;
  if (!cache_->is_unloading(&ran)) {
    ran = gen_uniform_();
    cache_->load(ran);
  }
  TRACE("ran: " << ran);
  return ran;
}

int Random::uniform(const int min, const int max) {
  ASSERT(max >= min, "max:" << max << " must be > min:" << min);
  return static_cast<int>(uniform() * (max - min + 1)) + min;
}

void Random::serialize_random_(std::ostream& ostr) const {
  feasst_serialize_version(979, ostr);
  feasst_serialize(cache_, ostr);
  feasst_serialize(is_seeded_, ostr);
}

std::map<std::string, std::shared_ptr<Random> >& Random::deserialize_map() {
  static std::map<std::string, std::shared_ptr<Random> >* ans =
     new std::map<std::string, std::shared_ptr<Random> >();
  return *ans;
}

std::shared_ptr<Random> Random::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

std::shared_ptr<Random> Random::factory(const std::string name,
                                        argtype * args) {
  return template_factory(deserialize_map(), name, args);
}

Random::Random(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 979, "mismatch version: " << version);
  // HWH for unknown reasons, the below does not work
  // feasst_deserialize(cache_, istr);
  {
    int existing;
    istr >> existing;
    if (existing != 0) {
      cache_ = std::make_shared<Cache>(istr);
    }
  }
  feasst_deserialize(&is_seeded_, istr);
}

bool Random::coin_flip() {
  if (uniform() > 0.5) {
    return true;
  }
  return false;
}

// acknowledge: http://stackoverflow.com/questions/440133/
// how-do-i-create-a-random-alpha-numeric-string-in-c
std::string Random::alpha_numeric(const int length) {
  std::string str;
  const std::string alphanum = "0123456789"
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
  for (int index = 0; index < length; ++index) {
    std::stringstream ss;
    ss << alphanum[uniform(0, alphanum.size() - 1)];
    str.append(ss.str());
  }
  return str;
}

double Random::uniform_real(const double min, const double max) {
  return (max - min)*uniform() + min;
}

void Random::position_in_cube(const int dimension,
    const double length,
    Position * position) {
  position->set_to_origin(dimension);
  for (int dim = 0; dim < dimension; ++dim) {
    position->set_coord(dim, (uniform() - 0.5)*length);
  }
}

Position Random::position_in_cube(const int dimension, const double length) {
  Position position;
  position_in_cube(dimension, length, &position);
  return position;
}

void Random::position_in_cuboid(const Position& side_length,
    Position * position) {
  if (position->size() == 0) position->set_to_origin(side_length.dimension());
  for (int dim = 0; dim < side_length.dimension(); ++dim) {
    position->set_coord(dim, (uniform() - 0.5)*side_length.coord(dim));
  }
}

Position Random::position_in_cuboid(const Position& side_length) {
  Position position;
  position_in_cuboid(side_length, &position);
  return position;
}

void Random::unit_sphere_surface(Position * position) {
  if (position->dimension() == 3) {
    // thanks to http://mathworld.wolfram.com/SpherePointPicking.html
    const double theta = 2*PI*uniform();
    const double phi = std::acos(2*uniform() - 1);
    position->set_from_spherical(1., theta, phi);
  } else if (position->dimension() == 2) {
    const double theta = 2*PI*uniform();
    position->set_from_spherical(1., theta);
  } else {
    ERROR("unrecognized dimension(" << position->dimension() << ")");
  }
}

int Random::index_from_cumulative_probability(
    const std::vector<double>& cumulative) {
  ASSERT(std::abs(cumulative.back() - 1) < 10.*NEAR_ZERO,
    "the cumulative distribution must end with a value of unity");
  if (static_cast<int>(cumulative.size()) == 1) return 0;
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
  INFO("This should never happen: " << feasst_str(cumulative) <<
       " ran " << random_uniform);
  // catch nearly impossible round off error
  if (std::abs(cumulative.back() - 1) <= NEAR_ZERO) {
    return static_cast<int>(cumulative.size()) - 1;
  }
  ERROR("This should never ever EVER happen!");
  return -1;
}

RotationMatrix Random::rotation(const int dimension, const double tunable) {
  Position axis;
  RotationMatrix rot_mat;
  rotation(dimension, &axis, &rot_mat, tunable);
  return rot_mat;
}

void Random::quaternion(Position * quat) {
  ASSERT(quat->size() == 4, "not implemented for size " << quat->size());
  double s1 = 2;
  double s2 = 2;
  while (s1 > 1) {
    quat->set_coord(0, uniform());
    quat->set_coord(0, 2*quat->coord(0) - 1);
    quat->set_coord(1, uniform());
    quat->set_coord(1, 2*quat->coord(1) - 1);
    s1 = quat->coord(0)*quat->coord(0) + quat->coord(1)*quat->coord(1);
  }
  while (s2 > 1) {
    quat->set_coord(2, uniform());
    quat->set_coord(2, 2*quat->coord(2) - 1);
    quat->set_coord(3, uniform());
    quat->set_coord(3, 2*quat->coord(3) - 1);
    s2 = quat->coord(2)*quat->coord(2) + quat->coord(3)*quat->coord(3);
  }
  s1 = std::sqrt((1 - s1)/s2);
  quat->set_coord(2, quat->coord(2)*s1);
  quat->set_coord(3, quat->coord(3)*s1);
}

void Random::rotation(const int dimension,
    Position * quaternion_or_axis,
    RotationMatrix * rot_mat,
    const double tunable) {
  const bool completely_random = std::abs(tunable + 1) < NEAR_ZERO;
  ASSERT(completely_random || tunable >= 0,
    "tunable: " << tunable << " must be -1 or > 0");
  if (rot_mat->num_rows() != dimension) rot_mat->set_size(dimension, dimension);
  if (dimension == 3) {
    if (quaternion_or_axis->dimension() != 3) {
      quaternion_or_axis->set_vector({0., 0., 0., 0.});
    }
    quaternion(quaternion_or_axis);
    if (!completely_random) {
      quaternion_or_axis->multiply(tunable);
      quaternion_or_axis->add_to_coord(3, 1.);
      const double size = quaternion_or_axis->distance();
      quaternion_or_axis->divide(size);
    }
    rot_mat->quaternion(*quaternion_or_axis);
  } else if (dimension == 2) {
    if (quaternion_or_axis->dimension() != 2) {
      quaternion_or_axis->set_vector({0., 0.});
    }
    double max_angle = tunable;
    if (completely_random) max_angle = 180.;
    const double angle = uniform_real(-max_angle, max_angle);
    rot_mat->axis_angle_opt(*quaternion_or_axis, angle);
  } else {
    ERROR("not implemented for dimension: " << dimension);
  }
}

void Random::position_in_spherical_shell(
    const double lower,
    const double upper,
    Position * position) {
  ASSERT(upper >= 0, "upper(" << upper << ") must be positive");
  ASSERT(lower >= 0, "lower(" << lower << ") must be positive");
  ASSERT(upper >= lower, "upper (" << upper << ") must be greater than "
         << "lower (" << lower << ")");
  const int dimen = position->dimension();
  if (dimen == 3) {
    const double lower3 = std::pow(lower, 3);
    const double upper3 = std::pow(upper, 3);
    const double factor = std::pow(uniform()*(upper3 - lower3) + lower3, 1./3.);
    unit_sphere_surface(position);
    position->multiply(factor);
  } else if (dimen == 2) {
    const double theta = 2*PI*uniform();
    const double lower_sq = lower*lower;
    const double upper_sq = upper*upper;
    const double factor = std::sqrt(uniform()*(upper_sq - lower_sq) + lower_sq);
    position->set_coord(0, factor*std::cos(theta));
    position->set_coord(1, factor*std::sin(theta));
  } else {
    FATAL("unrecognized dimension: " << dimen);
  }
}

double Random::standard_normal() {
  const double u = uniform();
  const double v = uniform();
  return std::sqrt(-2.*std::log(u))*std::cos(2.*PI*v);

// Frenkel and Smit a la Numerical Recipes
//  double v1, v2, r = 2.;
//  while (r >= 1.) {
//    v1 = 2.*uniform() - 1.;
//    v2 = 2.*uniform() - 1.;
//    r=v1*v1 + v2*v2;
//  }
//  return v1*std::sqrt(-2*std::log(r)/r);
}

// double Random::bond_angle(const double theta0,
//     const double beta_spring_constant,
//     const int exponent,
//     const int dimension,
//     const double minimum_angle) {
//   if (dimension == 2) {
//     FATAL("implmement flexible bonds in 2D.");
//   } else if (dimension != 3) {
//     FATAL("unrecognized dimension: " << dimension);
//   }
//   int attempt = 0;
//   while (attempt < 1e6) {
//     const double theta = minimum_angle + (PI - minimum_angle)*uniform();
//     const double dtheta = radians_to_degrees(theta - theta0);
//     const double beta_delta_U = beta_spring_constant*pow(dtheta, exponent);
//     if (uniform() < std::sin(theta)*std::exp(-beta_delta_U)) return theta;
//     ++attempt;
//   }
//   FATAL("max attempts reached");
// }

int Random::gen_uniform_(const int min, const int max) {
  return min + static_cast<int>(gen_uniform_()*(max - min));
}

const Cache& Random::cache() const { return *cache_; }

void Random::set_cache_to_load(const bool load) { cache_->set_load(load); }

void Random::set_cache_to_unload(const Random& random) {
  cache_->set_unload(random.cache());
}

}  // namespace feasst
