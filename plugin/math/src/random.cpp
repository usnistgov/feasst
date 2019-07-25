#include <string>
#include <sstream>
#include "math/include/random.h"
#include "utils/include/utils_io.h"

namespace feasst {

void seed_random_by_date() {
  const int t = time(NULL);
  srand ( t );
  INFO("time(seed): " << t);
}

void seed_random(const int seed) {
  srand ( seed );
  INFO("Initializing random number generator for reproduction with seed("
    << seed << ")");
}

int Random::uniform(const int min, const int max) {
  auto dis_int = std::uniform_int_distribution<int>(min, max);
  return dis_int(generator_);
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
  ASSERT(dimension == 3, "size error");
  position->set_to_origin_3D();
  for (int dim = 0; dim < dimension; ++dim) {
    position->set_coord(dim, (2.*uniform() - 1.)*length);
  }
}

Position Random::position_in_cube(const int dimension, const double length) {
  Position position;
  position_in_cube(dimension, length, &position);
  return position;
}

void Random::position_in_cuboid(const Position& side_length, Position * position) {
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

int Random::index_from_cumulative_probability(
    const std::vector<double>& cumulative) {
  ASSERT(std::abs(cumulative.back() - 1) < NEAR_ZERO,
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
  ERROR("should never reach this point");
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

void Random::reseed_() {
  const int seed = rand();
  TRACE("seed " << seed << " address " << this);
  generator_ = std::mt19937(seed);
  dis_double_ = std::uniform_real_distribution<double>(0.0, 1.0);
}

void Random::serialize(std::ostream& ostr) const {
  feasst_serialize_version(101, ostr);
  ostr << MAX_PRECISION;
  ostr << generator_;
}

Random::Random(std::istream& istr) {
  const int verison = feasst_deserialize_version(istr);
  ASSERT(verison == 101, "version");
  istr >> generator_;
}

}  // namespace feasst
