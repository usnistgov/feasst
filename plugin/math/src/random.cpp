#include <cmath>  // acos
#include <string>
#include <sstream>
#include "utils/include/utils_io.h"
#include "utils/include/serialize.h"
#include "math/include/random.h"
#include "math/include/constants.h"
#include "math/include/matrix.h"

namespace feasst {

Random::Random(const argtype& args) {}

// parsing seed in constructor leads to pure virtual function reseed_
void Random::parse_seed_(const argtype& args) {
  Arguments args_;
  args_.init(args);
  if (args_.key("seed").used()) {
    const std::string seed_str = args_.str();
    if (seed_str == "time") {
      seed_by_time();
    } else if (seed_str == "default") {
      seed();
    } else {
      args_.init(args);
      seed(args_.key("seed").integer());
    }
  }
}

void Random::seed_by_time() {
  const int t = time(NULL);
  srand(t);
  INFO("time(seed): " << t);
  reseed_();
  is_seeded_ = true;
}

void Random::seed(const int seed) {
  srand(seed);
  INFO("Initializing random number generator for reproduction with seed("
    << seed << ")");
  reseed_();
  is_seeded_ = true;
}

double Random::uniform() {
  if (!is_seeded_) {
    seed_by_time();
  }
  double ran;
  if (!cache_.is_unloading(&ran)) {
    ran = gen_uniform_();
    cache_.load(ran);
  }
  DEBUG("ran: " << ran);
  return ran;
}

int Random::uniform(const int min, const int max) {
  ASSERT(max >= min, "max:" << max << " must be > min:" << min);
  return static_cast<int>(uniform() * (max - min + 1)) + min;
}

void Random::serialize_random_(std::ostream& ostr) const {
  feasst_serialize_version(979, ostr);
  feasst_serialize_fstobj(cache_, ostr);
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

Random::Random(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 979, "mismatch version: " << version);
  feasst_deserialize_fstobj(&cache_, istr);
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
  ASSERT(side_length.dimension() == 3, "size error");
  position->set_to_origin_3D();
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
  ERROR("This should never happen even more!");
  return -1;
}

RotationMatrix Random::rotation(const int dimension, const double max_angle) {
  Position axis;
  if (dimension == 3) {
    axis.set_vector({0., 0., 0.});
  } else {
    ERROR("implement dimension: " << dimension);
  }
  unit_sphere_surface(&axis);
  const double angle = uniform_real(-max_angle, max_angle);
  return RotationMatrix().axis_angle(axis, angle);
}

}  // namespace feasst
