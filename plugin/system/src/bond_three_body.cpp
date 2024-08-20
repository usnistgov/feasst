
#include <cmath>
#include "utils/include/debug.h"
#include "utils/include/serialize_extra.h"
#include "math/include/utils_math.h"
#include "math/include/random.h"
#include "math/include/constants.h"
#include "math/include/position.h"
#include "configuration/include/bond.h"
#include "system/include/bond_three_body.h"

namespace feasst {

std::map<std::string, std::shared_ptr<BondThreeBody> >& BondThreeBody::deserialize_map() {
  static std::map<std::string, std::shared_ptr<BondThreeBody> >* ans =
     new std::map<std::string, std::shared_ptr<BondThreeBody> >();
  return *ans;
}

void BondThreeBody::serialize(std::ostream& ostr) const { FATAL("not implemented"); }

std::shared_ptr<BondThreeBody> BondThreeBody::create(std::istream& istr) const {
  FATAL("not implemented");
}

std::shared_ptr<BondThreeBody> BondThreeBody::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

void BondThreeBody::serialize_bond_three_body_(std::ostream& ostr) const {
  feasst_serialize_version(943, ostr);
}

BondThreeBody::BondThreeBody(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(943 == version, "mismatch version: " << version);
}

double BondThreeBody::energy(const Position& ri, const Position& rj,
    const Position& rk, const Bond& angle) const {
  return energy(rj.vertex_angle_radians(ri, rk), angle);
}

double BondThreeBody::random_angle_radians(const Angle& angle,
    const double beta, const int dimension, Random * random) const {
  ASSERT(dimension == 2 || dimension == 3,
    "unrecognized dimension: " << dimension);
  double min_rad = 0.;
  if (angle.has_property("minimum_degrees")) {
    min_rad = degrees_to_radians(angle.property("minimum_degrees"));
  }
  int attempt = 0;
  while (attempt < 1e6) {
    const double radians = min_rad + (PI - min_rad)*random->uniform();
    const double en = energy(radians, angle);
    double jacobian = 1.;
    if (dimension == 3) jacobian = std::sin(radians);
    if (random->uniform() < jacobian*std::exp(-beta*en)) {
      return radians;
    }
    ++attempt;
  }
  FATAL("max attempts reached");
}

void BondThreeBody::random_branch(
    const Angle& a2a1m1,
    const Angle& a2a1m2,
    const Angle& m1a1m2,
    const double beta,
    const bool is_position_held,
    double * radians_a2a1m1,
    double * radians_a2a1m2,
    double * radians_m1a1m2,
    Random * random,
    double * ln_met,
    const Position * const a1,
    const Position * const a2,
    const Position * const m1,
    const Position * const m2) const {
  if (is_position_held) return;
  ASSERT(a2a1m1.model() == a2a1m2.model() &&
         a2a1m1.model() == m1a1m2.model(), "Branch model mismatch: "
    << a2a1m1.model() << " " << a2a1m2.model() << " " << m1a1m2.model());
  const int dimen = 3;
  Position unit_m1(dimen), unit_m2(dimen);
  bool accept = false;
  int attempt = 0;
  while (!accept) {
    random->unit_sphere_surface(&unit_m1);
    random->unit_sphere_surface(&unit_m2);
    // assumes anchor 1 is origin and anchor 2 is unit vector along x-axis
    const double tm1 = std::acos(unit_m1.coord(0));
    const double tm2 = std::acos(unit_m2.coord(0));
    const double tm12 = std::acos(unit_m1.dot_product(unit_m2));
    if (random->uniform() < std::exp(-beta*(
        energy(tm1, a2a1m1) +
        energy(tm2, a2a1m2) +
        energy(tm12, m1a1m2)))) {
      accept = true;
      *radians_a2a1m1 = tm1;
      *radians_a2a1m2 = tm2;
      *radians_m1a1m2 = tm12;
    }
    ++attempt;
    if (attempt == 1e6) FATAL("max attempts reached");
  }
}

}  // namespace feasst
