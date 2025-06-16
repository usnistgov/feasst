#include <cmath>
#include "utils/include/debug.h"
#include "utils/include/io.h"
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/position.h"
#include "math/include/utils_math.h"
#include "math/include/constants.h"

namespace feasst {

double Position::coord(const int dimension) const {
  return coord_[dimension];
}

void Position::set_coord(const int dimension, const double coord) {
  ASSERT(dimension < size(), "dimension(" << dimension << ") < size(" << size()
    << ")");
  coord_[dimension] = coord;
}

void Position::add_to_coord(const int dimension, const double coord) {
  ASSERT(dimension < size(), "dimension(" << dimension << ") < size(" << size()
    << ")");
  coord_[dimension] += coord;
}

void Position::add(const Position &position) {
  ASSERT(position.size() == size(),
    "size:" << size() << " != position.size:" << position.size());
  for (int dim = 0; dim < size(); ++dim) {
    coord_[dim] += position.coord(dim);
  }
}

void Position::subtract(const Position &position) {
  ASSERT(position.size() == size(), "size mismatch");
  for (int dim = 0; dim < size(); ++dim) {
    coord_[dim] -= position.coord(dim);
  }
}

void Position::divide(const Position &position) {
  ASSERT(position.size() == size(), "size mismatch");
  for (int dim = 0; dim < size(); ++dim) {
    coord_[dim] /= position.coord(dim);
  }
}

double Position::dot_product(const Position &position) const {
  return dot_product(position.coord());
}

double Position::dot_product(const std::vector<double>& vec) const {
  ASSERT(static_cast<int>(vec.size()) == size(), "size mismatch");
  double prod = 0.;
  for (int dim = 0; dim < size(); ++dim) {
    prod += coord_[dim]*vec[dim];
  }
  return prod;
}

double Position::squared_distance() const {
  double dist_sq = 0;
  for (const double& coord : coord_) {
    dist_sq += coord*coord;
  }
  return dist_sq;
}

double Position::distance() const {
  return std::sqrt(squared_distance());
}

double Position::squared_distance(const Position& position) const {
  ASSERT(dimension() == position.dimension(), "Err");
  double dist = 0;
  for (int dim = 0; dim < dimension(); ++dim) {
    const double diff = coord_[dim] - position.coord_[dim];
    dist += diff*diff;
  }
  return dist;
}

double Position::distance(const Position& position) const {
  return std::sqrt(squared_distance(position));
}

std::string Position::str() const {
  return feasst_str(coord_, true);
}

void Position::set_to_origin() {
  for (double& value : coord_) {
    value = 0.;
  }
}

void Position::set_to_origin(const int dimension) {
  coord_.resize(dimension);
  set_to_origin();
}

double Position::cosine(const Position& position) const {
  const double dist = position.distance();
  ASSERT(dist > NEAR_ZERO, "dist: " << dist << " is too small. "
    << "May be caused by particles in an angle on top of each other.");
  double cos = dot_product(position)/distance()/dist;
  ASSERT(std::abs(cos) < 1 + NEAR_ZERO, "|cos: " << cos << "| > 1");
  if (cos < -1) cos = -1.;
  if (cos > 1) cos = 1.;
  return cos;
}

void Position::divide(const double denominator) {
  for (double& coord : coord_) {
    coord /= denominator;
  }
}

void Position::normalize() {
  const double dist = distance();
  ASSERT(std::abs(dist) > NEAR_ZERO, " cannot normalize a 0 vector: " << str());
  divide(dist);
}

bool Position::is_equal(const Position& position,
                        const double tolerance) const {
  if (size() != position.size()) {
    return false;
  }
  for (int index = 0; index < static_cast<int>(size()); ++index) {
    if (std::abs(position.coord(index) - coord_[index]) > tolerance) {
      return false;
    }
  }
  return true;
}

bool Position::is_equal(const Position& position) const {
  return is_equal(position, NEAR_ZERO);
}

Position& Position::set_from_spherical(const std::vector<double> &vec) {
  const int dimen = static_cast<int>(vec.size());
  if (dimen == 3) {
    set_from_spherical(vec[0], vec[1], vec[2]);
  } else if (dimen == 2) {
    set_from_spherical(vec[0], vec[1]);
  } else {
    ERROR("unrecognized dimension:" << vec.size());
  }
  return *this;
}

void Position::set_from_spherical(const double rho,
                                  const double theta,
                                  const double phi) {
  const double sine_phi = sin(phi);
  set_vector({
    rho*sine_phi*cos(theta),
    rho*sine_phi*sin(theta),
    rho*cos(phi)});
}

void Position::set_from_spherical(const double rho,
                                  const double theta) {
  set_vector({rho*cos(theta), rho*sin(theta)});
}

void Position::multiply(const double constant) {
  for (double& coord : coord_) {
    coord *= constant;
  }
}

void Position::serialize(std::ostream& sstr) const {
  feasst_serialize_version(3914, sstr);
  feasst_serialize(coord_, sstr);
}

Position::Position(std::istream& sstr) {
  const int version = feasst_deserialize_version(sstr);
  ASSERT(version == 3914, "unrecognized verison: " << version);
  feasst_deserialize(&coord_, sstr);
}

Position Position::cross_product(const Position& position) const {
  ASSERT(dimension() == 3 && position.dimension() == 3,
    "implemented for 3D only.");
  return Position().set_vector({
    coord(1) * position.coord(2) - coord(2) * position.coord(1),
    coord(2) * position.coord(0) - coord(0) * position.coord(2),
    coord(0) * position.coord(1) - coord(1) * position.coord(0)});
}

double Position::nearest_distance_to_axis(const Position& point1,
                                          const Position& point2) const {
  ASSERT(dimension() == 3 && point1.dimension() == 3 && point2.dimension() == 3,
    "implemented for 3D only.");
  const double p00 = coord_[0],
    p01 = coord_[1],
    p02 = coord_[2],
    p10 = point1.coord_[0],
    p11 = point1.coord_[1],
    p12 = point1.coord_[2],
    p20 = point2.coord_[0],
    p21 = point2.coord_[1],
    p22 = point2.coord_[2],
    // 3 relative vectors
    d010 = p00 - p10,
    d011 = p01 - p11,
    d012 = p02 - p12,
    d020 = p00 - p20,
    d021 = p01 - p21,
    d022 = p02 - p22,
    d210 = p20 - p10,
    d211 = p21 - p11,
    d212 = p22 - p12,
    // cross product
    c0 = d011*d022 - d012*d021,
    c1 = d012*d020 - d010*d022,
    c2 = d010*d021 - d011*d020;
  return std::sqrt((c0*c0 + c1*c1 + c2*c2)/(d210*d210 + d211*d211 + d212*d212));

//  // This is the old, non optimized version
//  Position p01 = *this, p02 = *this, p21 = point2;
//  p01.subtract(point1);
//  p02.subtract(point2);
//  p21.subtract(point1);
//  return ((p01.cross_product(p02)).distance())/p21.distance();
}

Position::Position(argtype * args) {
  double x = 0, y = 0, z = 0;
  if (!used("x", *args)) {
    const std::string csv = feasst::str("csv", args, "");
    if (!csv.empty()) {
      std::vector<double> vals;
      for (const std::string& x : split(csv, ',')) {
        vals.push_back(str_to_double(x));
      }
      set_vector(vals);
    } else {
      const int dimension = integer("dimension", args, 0);
      ASSERT(dimension >= 0, "dimension cannot be negative.");
      coord_.resize(dimension, 0.);
    }
  } else {
    WARN("Deprecated Position::x->csv");
    x = dble("x", args);
    if (used("y", *args)) {
      y = dble("y", args);
      if (used("z", *args)) {
        z = dble("z", args);
        set_vector({x, y, z});
      } else {
        set_vector({x, y});
      }
    } else {
      set_vector({x});
    }
  }
}
Position::Position(argtype args) : Position(&args) {
  feasst_check_all_used(args);
}

void Position::orthogonal(const Position& orthogonal) {
  ASSERT(orthogonal.size() == 3,
    "only implemented for 3D");
  if (size() != 3) set_to_origin(3);

  // find index with least absolute value
  int min = 0;
  if (std::abs(orthogonal.coord(1)) < std::abs(orthogonal.coord(min))) min = 1;
  if (std::abs(orthogonal.coord(2)) < std::abs(orthogonal.coord(min))) min = 2;

  // define orthogonal vector by setting the min to 0 and swapping the other
  // two indices, setting one of the swapped negative.
  if (min == 0) {
    coord_[0] = 0.;
    coord_[1] = -orthogonal.coord(2);
    coord_[2] = orthogonal.coord(1);
  } else if (min == 1) {
    coord_[0] = -orthogonal.coord(2);
    coord_[1] = 0.;
    coord_[2] = orthogonal.coord(0);
  } else {
    coord_[0] = -orthogonal.coord(1);
    coord_[1] = orthogonal.coord(0);
    coord_[2] = 0.;
  }
  normalize();
  ASSERT(std::abs(dot_product(orthogonal)) < 100*NEAR_ZERO,
    "self(" << str() << " is not orthogonal with: " << orthogonal.str() <<
    " if the dot product of the two is: " << dot_product(orthogonal));
}

Position& Position::set_vector(const std::vector<double> &vec) {
  coord_ = vec;
  return *this;
}

void Position::reflect(const Position& reflection_point) {
  for (int dim = 0; dim < dimension(); ++dim) {
    coord_[dim] = 2.*reflection_point.coord(dim) - coord_[dim];
  }
}

Position Position::spherical() const {
  Position sph = *this;
  spherical(&sph);
  return sph;
}

void Position::spherical(Position * result) const {
  ASSERT(dimension() == 2 || dimension() == 3,
    "unrecognized dimension: " << dimension());
  const double r = distance();
  result->set_coord(0, r);
  if (r > NEAR_ZERO) {
    result->set_coord(1, std::atan2(coord(1), coord(0)));
  } else {
    result->set_coord(1, 0.);
  }
  if (dimension() == 3) {
    if (r > NEAR_ZERO) {
      result->set_coord(2, std::acos(coord(2)/r));
    } else {
      result->set_coord(2, 0.);
    }
  }
}

double Position::vertex_angle_radians(const Position& ri,
                                      const Position& rk) const {
  Position rij = ri;
  rij.subtract(*this);
  Position rkj = rk;
  rkj.subtract(*this);
  double rad = std::acos(rij.cosine(rkj));
  if (ri.dimension() == 2) {
    // if the z-dimension of the cross product is positive, reverse
    if (rij.coord(0)*rkj.coord(1) > rij.coord(1)*rkj.coord(0)) {
      rad = 2.*PI - rad;
    }
  }
  return rad;
}

double Position::torsion_angle_radians(const Position& rj, const Position& rk,
    const Position& rl) const {
  Position rij = *this;
  rij.subtract(rj);
  Position rjk = rj;
  rjk.subtract(rk);
  Position rkl = rk;
  rkl.subtract(rl);
  Position n1 = rkl.cross_product(rjk);
  const double n1_mag = n1.distance();
  ASSERT(std::abs(n1_mag) > NEAR_ZERO,
    "n1: " << n1.str() << " is too small. " << str() << " " << rj.str() <<
    " " << rk.str() << " " << rl.str());
  DEBUG("n1 " << n1.str());
  Position n2 = rjk.cross_product(rij);
  const double n2_mag = n2.distance();
  ASSERT(std::abs(n2_mag) > NEAR_ZERO, "n2 is too small");
  DEBUG("n2 " << n2.str());
  double cos = n1.dot_product(n2)/n1_mag/n2_mag;
  if (cos < -1) cos = -1.;
  if (cos > 1) cos = 1.;
  const double angle = std::acos(cos);
  ASSERT(!std::isnan(angle), "angle: " << angle << " cos " << cos << " n1 "
    << n1.str() << " n2 " << n2.str());
  return angle;
}

}  // namespace feasst
