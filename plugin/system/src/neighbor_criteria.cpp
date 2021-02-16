#include <cmath>
#include "system/include/neighbor_criteria.h"
#include "utils/include/serialize.h"
#include "utils/include/debug.h"
#include "math/include/utils_math.h"
#include "math/include/constants.h"
#include "math/include/random.h"
#include "configuration/include/domain.h"

namespace feasst {

NeighborCriteria::NeighborCriteria(argtype args) {
  reference_potential_ = integer("reference_potential", &args, -1);
  potential_index_ = integer("potential_index", &args, 0);
  energy_maximum_ = dble("energy_maximum", &args, -NEAR_ZERO);
  DEBUG("energy_maximum " << energy_maximum_);
  ASSERT(energy_maximum_ < 0.,
    "energy_maximum:" << energy_maximum_ << " must be less than zero. " <<
    "Otherwise, self interactions and particles outside of cutoff will be added");
  minimum_distance_sq_ = std::pow(dble("minimum_distance", &args, 0), 2);
  maximum_distance_sq_ = std::pow(
    dble("maximum_distance", &args, std::sqrt(NEAR_INFINITY)), 2);
  site_type0_ = integer("site_type0", &args, -1);
  site_type1_ = integer("site_type1", &args, -1);
  ASSERT((site_type0_ == -1 && site_type1_ == -1) ||
         (site_type0_ != -1 && site_type1_ != -1),
    "site_type0: " << site_type0_ << " and site_type1: " << site_type1_ <<
    " must be either both -1 or neither -1");
  check_all_used(args);
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

double NeighborCriteria::volume(const int dimension) const {
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
  DEBUG("squared_distance " << squared_distance);
  return squared_distance > minimum_distance_sq_ &&
         squared_distance < maximum_distance_sq_;
}

}  // namespace feasst
