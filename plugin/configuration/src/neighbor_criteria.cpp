#include <cmath>
#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "utils/include/debug.h"
#include "math/include/utils_math.h"
#include "math/include/constants.h"
#include "math/include/random.h"
#include "configuration/include/domain.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/neighbor_criteria.h"

namespace feasst {

NeighborCriteria::NeighborCriteria(argtype * args) {
  if (used("reference_potential", *args)) {
    FATAL("Deprecated NeighborCriteria::reference_potential->ref.");
  }
  reference_potential_ = integer("reference_potential", args, -1);
  ref_ = str("ref", args, "");
  potential_index_ = integer("potential_index", args, 0);
  energy_maximum_ = dble("energy_maximum", args,
                         std::numeric_limits<double>::max());
  DEBUG("energy_maximum " << energy_maximum_);
  minimum_distance_sq_ = std::pow(dble("minimum_distance", args, 0), 2);
  maximum_distance_sq_ = std::pow(
    dble("maximum_distance", args, std::sqrt(NEAR_INFINITY)), 2);
  site_type0_name_ = str("site_type0", args, "-1");
  site_type1_name_ = str("site_type1", args, "-1");
  ASSERT((site_type0_name_ == "-1" && site_type1_name_ == "-1") ||
         (site_type0_name_ != "-1" && site_type1_name_ != "-1"),
    "site_type0: " << site_type0_name_ << " and site_type1: " << site_type1_name_ <<
    " must be either both -1 or neither -1");
  site_type0_alt_name_ = str("site_type0_alt", args, "-1");
  site_type1_alt_name_ = str("site_type1_alt", args, "-1");
  ASSERT((site_type0_alt_name_ == "-1" && site_type1_alt_name_ == "-1") ||
         (site_type0_alt_name_ != "-1" && site_type1_alt_name_ != "-1"),
    "site_type0_alt: " << site_type0_alt_name_ << " and site_type1_alt: " << site_type1_alt_name_ <<
    " must be either both -1 or neither -1");

  // Make a best guess for integer site types
  try {
    site_type0_ = std::stoi(site_type0_name_);
  } catch (...) {
    site_type0_ = -1;
  }
  try {
    site_type1_ = std::stoi(site_type1_name_);
  } catch (...) {
    site_type1_ = -1;
  }
  try {
    site_type0_alt_ = std::stoi(site_type0_alt_name_);
  } catch (...) {
    site_type0_alt_ = -1;
  }
  try {
    site_type1_alt_ = std::stoi(site_type1_alt_name_);
  } catch (...) {
    site_type1_alt_ = -1;
  }
}
NeighborCriteria::NeighborCriteria(argtype args) : NeighborCriteria(&args) {
  feasst_check_all_used(args);
}

void NeighborCriteria::name_to_index(const ParticleFactory& unique_types) {
  DEBUG("site_type0_name_ " << site_type0_name_);
  DEBUG("site_type1_name_ " << site_type1_name_);
  DEBUG("site_type0_alt_name_ " << site_type0_alt_name_);
  DEBUG("site_type1_alt_name_ " << site_type1_alt_name_);
  if (site_type0_name_ == "-1" && site_type1_name_ == "-1") {
    site_type0_ = site_type1_ = -1;
  } else if (site_type0_name_ != "-1" && site_type1_name_ != "-1") {
    site_type0_ = unique_types.site_type_name_to_index(site_type0_name_);
    site_type1_ = unique_types.site_type_name_to_index(site_type1_name_);
  } else {
    FATAL("site_type0_name_:" << site_type0_name_ << " and "
      << "site_type1_name_:" << site_type1_name_ << " should either "
      << "both be -1 or neither");
  }
  if (site_type0_alt_name_ == "-1" && site_type1_alt_name_ == "-1") {
    site_type0_alt_ = site_type1_alt_ = -1;
  } else if (site_type0_alt_name_ != "-1" && site_type1_alt_name_ != "-1") {
    site_type0_alt_ = unique_types.site_type_name_to_index(site_type0_alt_name_);
    site_type1_alt_ = unique_types.site_type_name_to_index(site_type1_alt_name_);
  } else {
    FATAL("site_type0_alt_name_:" << site_type0_alt_name_ << " and "
      << "site_type1_alt_name_:" << site_type1_alt_name_ << " should either "
      << "both be -1 or neither");
  }
}

bool NeighborCriteria::is_accepted(const double energy,
                                   const double squared_distance,
                                   const int site_type0,
                                   const int site_type1) const {
  const bool is_distance = squared_distance > minimum_distance_sq_ &&
                           squared_distance < maximum_distance_sq_;
  DEBUG("is_distance " << is_distance);
  if (!is_distance) {
    return false;
  }
  const bool is_energy = energy < energy_maximum_;
  DEBUG("is_energy " << is_energy << " energy: " << energy);
  if (!is_energy) {
    return false;
  }
  const bool is_type = site_type0_ == -1 ||
    ( (site_type0_ == site_type0 && site_type1_ == site_type1) ||
      (site_type0_ == site_type1 && site_type1_ == site_type0) ||
      (site_type0_alt_ != -1 &&
        ( (site_type0_alt_ == site_type0 && site_type1_alt_ == site_type1) ||
          (site_type0_alt_ == site_type1 && site_type1_alt_ == site_type0) ) ));
  DEBUG("is_type " << is_type << " site_type0 " << site_type0 <<
    " site_type0_ " << site_type0_ << " site_type1 " << site_type1 <<
    " site_type1_ " << site_type1_ << " site_type0_alt_ " << site_type0_alt_ <<
    " site_type1_alt_ " << site_type1_alt_);
  if (!is_type) {
    return false;
  }
  return true;
}

void NeighborCriteria::serialize(std::ostream& ostr) const {
  feasst_serialize_version(906, ostr);
  feasst_serialize(reference_potential_, ostr);
  feasst_serialize(ref_, ostr);
  feasst_serialize(potential_index_, ostr);
  feasst_serialize(energy_maximum_, ostr);
  feasst_serialize(minimum_distance_sq_, ostr);
  feasst_serialize(maximum_distance_sq_, ostr);
  feasst_serialize(site_type0_, ostr);
  feasst_serialize(site_type1_, ostr);
  feasst_serialize(site_type0_alt_, ostr);
  feasst_serialize(site_type1_alt_, ostr);
  feasst_serialize(site_type0_name_, ostr);
  feasst_serialize(site_type1_name_, ostr);
  feasst_serialize(site_type0_alt_name_, ostr);
  feasst_serialize(site_type1_alt_name_, ostr);
}

NeighborCriteria::NeighborCriteria(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 903 && version <= 906, "version: " << version);
  feasst_deserialize(&reference_potential_, istr);
  if (version >= 906) {
    feasst_deserialize(&ref_, istr);
  }
  feasst_deserialize(&potential_index_, istr);
  feasst_deserialize(&energy_maximum_, istr);
  feasst_deserialize(&minimum_distance_sq_, istr);
  feasst_deserialize(&maximum_distance_sq_, istr);
  feasst_deserialize(&site_type0_, istr);
  feasst_deserialize(&site_type1_, istr);
  if (version >= 904) {
    feasst_deserialize(&site_type0_alt_, istr);
    feasst_deserialize(&site_type1_alt_, istr);
  }
  if (version >= 905) {
    feasst_deserialize(&site_type0_name_, istr);
    feasst_deserialize(&site_type1_name_, istr);
    feasst_deserialize(&site_type0_alt_name_, istr);
    feasst_deserialize(&site_type1_alt_name_, istr);
  }
}

double NeighborCriteria::minimum_distance() const {
  return std::sqrt(minimum_distance_sq_);
}

double NeighborCriteria::maximum_distance() const {
  return std::sqrt(maximum_distance_sq_);
}

double NeighborCriteria::volume(const int dimension) const {
  const double volume = spherical_shell_volume(std::sqrt(minimum_distance_sq_),
    std::sqrt(maximum_distance_sq_),
    dimension);
  ASSERT(volume < NEAR_INFINITY, "NeighborCriteria volume is infinite. "
    << "A maximum_distance argument is required to compute volume.");
  return volume;
}

bool NeighborCriteria::is_position_accepted(
    const Position& position,
    const Domain& domain) {
  double squared_distance;
  if (!origin_) {
    origin_ = std::make_shared<Position>();
    rel_ = std::make_shared<Position>();
    pbc_ = std::make_shared<Position>();
  }
  if (origin_->dimension() == 0) {
    origin_->set_to_origin(domain.dimension());
    rel_->set_to_origin(domain.dimension());
    pbc_->set_to_origin(domain.dimension());
  }
  domain.wrap_opt(position, *origin_, rel_.get(), pbc_.get(),
                  &squared_distance);
  DEBUG("squared_distance " << squared_distance);
  return squared_distance > minimum_distance_sq_ &&
         squared_distance < maximum_distance_sq_;
}

int NeighborCriteria::reference_potential() const {
  //FATAL("used here");
  return reference_potential_;
}

}  // namespace feasst
