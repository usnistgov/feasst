
#include "system/include/system.h"
#include "utils/include/debug.h"

namespace feasst {

void System::add(const Configuration& configuration) {
  configurations_.push_back(configuration);
}

const Configuration& System::configuration(const int config) const {
  return configurations_[config];
}

Configuration* System::get_configuration(const int config) {
  return &configurations_[config];
}

int System::dimension(const int config) const {
  const int dim = configurations_[0].dimension();
  for (const Configuration& config : configurations_) {
    ASSERT(dim == config.dimension(), "dimensions of configs do not match");
  }
  return dim;
}

void System::add_to_unoptimized(const Potential& potential) {
  // HWH assume one config
  unoptimized_.add_potential(potential);
  unoptimized_.precompute(unoptimized_.num() - 1, &configurations_[0]);
}

void System::set_unoptimized(const int index, const Potential& potential) {
  // HWH assume one config
  unoptimized_.set_potential(index, potential);
  unoptimized_.precompute(unoptimized_.num() - 1, &configurations_[0]);
}

void System::add_to_optimized(const Potential& potential) {
  is_optimized_ = true;
  // HWH assume one config
  optimized_.add_potential(potential);
  optimized_.precompute(optimized_.num() - 1, &configurations_[0]);
}

PotentialFactory * System::reference_(const int index) {
  ASSERT(index < static_cast<int>(references_.size()),
    "unrecognized reference");
  return &references_[index];
}

void System::add_to_reference(const Potential& ref, const int index) {
  if (index == 0 and references_.size() == 0) {
    references_.push_back(PotentialFactory());
  }
  // HWH assume one config
  reference_(index)->add_potential(ref);
  reference_(index)->precompute(reference_(index)->num() - 1,
                                &configurations_[0]);
}

const Potential& System::reference(const int ref,
    const int potential) const {
  return references_[ref].potentials()[potential];
}

void System::precompute() {
  unoptimized_.precompute(&configurations_.front());
  if (is_optimized_) {
    optimized_.precompute(&configurations_.front());
  }
  for (PotentialFactory& ref : references_) {
    ref.precompute(&configurations_.front());
  }
}

double System::unoptimized_energy(const int config) {
  return unoptimized_.energy(&configurations_[0]);
}

PotentialFactory * System::potentials_() {
  if (is_optimized_) {
    return &optimized_;
  }
  return &unoptimized_;
}

double System::energy(const int config) {
  return potentials_()->energy(&configurations_[0]);
}

double System::energy(const Select& select, const int config) {
  return potentials_()->energy(select, &configurations_[0]);
}

double System::reference_energy(const Select& select, const int ref,
  const int config) {
  return reference_(ref)->energy(select, &configurations_[0]);
}

void System::serialize(std::ostream& sstr) const {
  feasst_serialize_version(1, sstr);
  feasst_serialize_fstobj(configurations_, sstr);
  unoptimized_.serialize(sstr);
  optimized_.serialize(sstr);
  feasst_serialize(is_optimized_, sstr);
  feasst_serialize_fstobj(references_, sstr);
}

System::System(std::istream& sstr) {
  feasst_deserialize_version(sstr);
  feasst_deserialize_fstobj(&configurations_, sstr);
  unoptimized_ = PotentialFactory(sstr);
  optimized_ = PotentialFactory(sstr);
  feasst_deserialize(&is_optimized_, sstr);
  feasst_deserialize_fstobj(&references_, sstr);
}

const PotentialFactory * System::const_potentials_() const {
  if (is_optimized_) {
    return &optimized_;
  }
  return &unoptimized_;
}

void System::load_cache(const bool load) {
  unoptimized_.load_cache(load);
  optimized_.load_cache(load);
  for (PotentialFactory& ref : references_) {
    ref.load_cache(load);
  }
}

void System::unload_cache(const System& system) {
  unoptimized_.unload_cache(system.unoptimized_);
  optimized_.unload_cache(system.optimized_);
  ASSERT(references_.size() == system.references_.size(), "size mismatch");
  for (int iref = 0; iref < static_cast<int>(references_.size()); ++iref) {
    references_[iref].unload_cache(system.references_[iref]);
  }
}

}  // namespace feasst
