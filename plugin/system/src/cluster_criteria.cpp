#include <cmath>
#include "system/include/cluster_criteria.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"

namespace feasst {

ClusterCriteria::ClusterCriteria(const argtype& args) {
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
}

bool ClusterCriteria::is_accepted(const double energy,
                                  const double squared_distance) const {
  if (energy < energy_maximum_ &&
      squared_distance > minimum_distance_sq_ &&
      squared_distance < maximum_distance_sq_) {
    return true;
  }
  return false;
}

void ClusterCriteria::serialize(std::ostream& ostr) const {
  feasst_serialize_version(903, ostr);
  feasst_serialize(reference_potential_, ostr);
  feasst_serialize(potential_index_, ostr);
  feasst_serialize(energy_maximum_, ostr);
  feasst_serialize(minimum_distance_sq_, ostr);
  feasst_serialize(maximum_distance_sq_, ostr);
}

ClusterCriteria::ClusterCriteria(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 903, "version");
  feasst_deserialize(&reference_potential_, istr);
  feasst_deserialize(&potential_index_, istr);
  feasst_deserialize(&energy_maximum_, istr);
  feasst_deserialize(&minimum_distance_sq_, istr);
  feasst_deserialize(&maximum_distance_sq_, istr);
}

double ClusterCriteria::minimum_distance() const {
  return std::sqrt(minimum_distance_sq_);
}

double ClusterCriteria::maximum_distance() const {
  return std::sqrt(maximum_distance_sq_);
}

}  // namespace feasst
