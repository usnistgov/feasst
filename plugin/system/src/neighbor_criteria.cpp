#include <cmath>
#include "system/include/neighbor_criteria.h"
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "math/include/utils_math.h"
#include "math/include/constants.h"
#include "math/include/random.h"
#include "configuration/include/domain.h"

namespace feasst {

NeighborCriteria::NeighborCriteria(const argtype& args) {
  Arguments args_(args);
  reference_potential_ = args_.key("reference_potential").dflt("-1").integer();
  potential_index_ = args_.key("potential_index").dflt("0").integer();
  energy_maximum_ = args_.key("energy_maximum").dflt(str(-NEAR_ZERO)).dble();
  DEBUG("energy_maximum " << energy_maximum_);
  ASSERT(energy_maximum_ < 0.,
    "energy_maximum:" << energy_maximum_ << " must be less than zero. " <<
    "Otherwise, self interactions and particles outside of cutoff will be added");
  minimum_distance_sq_ = std::pow(args_.key("minimum_distance").dflt("0").dble(), 2);
  maximum_distance_sq_ = std::pow(
    args_.key("maximum_distance").dflt(str(std::sqrt(NEAR_INFINITY))).dble(), 2);
  site_type0_ = args_.key("site_type0").dflt("-1").integer();
  site_type1_ = args_.key("site_type1").dflt("-1").integer();
  ASSERT((site_type0_ == -1 && site_type1_ == -1) ||
         (site_type0_ != -1 && site_type1_ != -1),
    "site_type0: " << site_type0_ << " and site_type1: " << site_type1_ <<
    " must be either both -1 or neither -1");
}

bool NeighborCriteria::is_accepted(const double energy,
                                   const double squared_distance,
                                   const int site_type0,
                                   const int site_type1) const {
  const bool is_distance = squared_distance > minimum_distance_sq_ &&
                           squared_distance < maximum_distance_sq_;
  const bool is_energy = energy < energy_maximum_;
  const bool is_type = site_type0_ == -1 ||
   ( (site_type0_ == site_type0 && site_type1_ == site_type1) ||
     (site_type0_ == site_type1 && site_type1_ == site_type0));
  if (is_distance && is_energy && is_type) {
    return true;
  }
  return false;
}

void NeighborCriteria::serialize(std::ostream& ostr) const {
  feasst_serialize_version(903, ostr);
  feasst_serialize(reference_potential_, ostr);
  feasst_serialize(potential_index_, ostr);
  feasst_serialize(energy_maximum_, ostr);
  feasst_serialize(minimum_distance_sq_, ostr);
  feasst_serialize(maximum_distance_sq_, ostr);
  feasst_serialize(site_type0_, ostr);
  feasst_serialize(site_type1_, ostr);
}

NeighborCriteria::NeighborCriteria(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 903, "version");
  feasst_deserialize(&reference_potential_, istr);
  feasst_deserialize(&potential_index_, istr);
  feasst_deserialize(&energy_maximum_, istr);
  feasst_deserialize(&minimum_distance_sq_, istr);
  feasst_deserialize(&maximum_distance_sq_, istr);
  feasst_deserialize(&site_type0_, istr);
  feasst_deserialize(&site_type1_, istr);
}

double NeighborCriteria::minimum_distance() const {
  return std::sqrt(minimum_distance_sq_);
}

double NeighborCriteria::maximum_distance() const {
  return std::sqrt(maximum_distance_sq_);
}

double NeighborCriteria::volume(const int dimension) {
  return spherical_shell_volume(std::sqrt(minimum_distance_sq_),
    std::sqrt(maximum_distance_sq_),
    dimension);
}

bool NeighborCriteria::is_position_accepted(
    const Position& position,
    const Domain& domain) {
  double squared_distance;
  if (origin_.dimension() == 0) {
    origin_.set_to_origin(domain.dimension());
    rel_.set_to_origin(domain.dimension());
    pbc_.set_to_origin(domain.dimension());
  }
  domain.wrap_opt(position, origin_, &rel_, &pbc_, &squared_distance);
  return squared_distance > minimum_distance_sq_ &&
         squared_distance < maximum_distance_sq_;
}

}  // namespace feasst
