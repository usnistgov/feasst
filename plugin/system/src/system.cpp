#include <cmath>  // isinf and isnan
#include "system/include/system.h"
#include "utils/include/debug.h"
#include "utils/include/serialize.h"

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
  unoptimized_.add(potential);
  unoptimized_.precompute(unoptimized_.num() - 1, &configurations_[0]);
}

void System::set_unoptimized(const int index, const Potential& potential) {
  // HWH assume one config
  unoptimized_.set(index, potential);
  unoptimized_.precompute(unoptimized_.num() - 1, &configurations_[0]);
}

void System::add_to_optimized(const Potential& potential) {
  is_optimized_ = true;
  // HWH assume one config
  optimized_.add(potential);
  optimized_.precompute(optimized_.num() - 1, &configurations_[0]);
}

PotentialFactory * System::reference_(const int index) {
  ASSERT(index < static_cast<int>(references_.size()),
    "unrecognized reference: " << index <<
    ". There are " << references_.size());
  return &references_[index];
}

void System::add_to_reference(const Potential& ref, const int index) {
  if (index == 0 and references_.size() == 0) {
    references_.push_back(PotentialFactory());
  }
  // HWH assume one config
  reference_(index)->add(ref);
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
  const double en = unoptimized_.energy(&configurations_[config]);
  ASSERT(!std::isinf(en) && !std::isnan(en),
    "Energy(" << en << ") is infinite or not "
    << "a number. Are particles on top of each other?");
  unoptimized_.finalize(configurations_[config].selection_of_all());
  ref_used_last_ = -1;
  DEBUG("ref_used_last_ " << ref_used_last_);
  return en;
}

PotentialFactory * System::potentials_() {
  if (is_optimized_) {
    return &optimized_;
  }
  return &unoptimized_;
}

double System::energy(const int config) {
  const double en = potentials_()->energy(&configurations_[config]);
  finalize(config);
  ref_used_last_ = -1;
  DEBUG("ref_used_last_ " << ref_used_last_);
  return en;
}

double System::perturbed_energy(const Select& select, const int config) {
  ref_used_last_ = -1;
  DEBUG("ref_used_last_ " << ref_used_last_);
  return potentials_()->energy(select, &configurations_[config]);
}

double System::reference_energy(const int ref, const int config) {
  ref_used_last_ = ref;
  DEBUG("ref_used_last_ " << ref_used_last_);
  return reference_(ref)->energy(&configurations_[config]);
}

double System::reference_energy(const Select& select,
    const int ref,
    const int config) {
  ref_used_last_ = ref;
  DEBUG("ref_used_last_ " << ref_used_last_);
  return reference_(ref)->energy(select, &configurations_[config]);
}

void System::serialize(std::ostream& sstr) const {
  feasst_serialize_version(7349, sstr);
  feasst_serialize_fstobj(configurations_, sstr);
  unoptimized_.serialize(sstr);
  optimized_.serialize(sstr);
  feasst_serialize(is_optimized_, sstr);
  feasst_serialize_fstobj(references_, sstr);
  feasst_serialize_endcap("System", sstr);
}

System::System(std::istream& sstr) {
  const int version = feasst_deserialize_version(sstr);
  ASSERT(version == 7349, "unrecognized verison: " << version);
  feasst_deserialize_fstobj(&configurations_, sstr);
  unoptimized_ = PotentialFactory(sstr);
  optimized_ = PotentialFactory(sstr);
  feasst_deserialize(&is_optimized_, sstr);
  feasst_deserialize_fstobj(&references_, sstr);
  feasst_deserialize_endcap("System", sstr);
}

const PotentialFactory& System::potentials() const {
  if (is_optimized_) {
    return optimized_;
  }
  return unoptimized_;
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

void System::finalize(const Select& select, const int config) {
  if (select.trial_state() == 2) {
    // finalize removal
    configurations_[config].remove_particles(select);
  }
  DEBUG("ref_used_last_ " << ref_used_last_);
  if (ref_used_last_ != -1) {
    references_[ref_used_last_].finalize(select);
  } else {
    potentials_()->finalize(select);
  }
}

void System::revert(const Select& select, const int config) {
  if (select.trial_state() == 3) {
    // revert addition
    configurations_[config].remove_particles(select);
  }
  DEBUG("ref_used_last_ " << ref_used_last_);
  if (ref_used_last_ != -1) {
    references_[ref_used_last_].revert(select);
  } else {
    potentials_()->revert(select);
  }
}

void System::check() const {
  unoptimized_.check();
  optimized_.check();
  for (const PotentialFactory& ref : references_) {
    ref.check();
  }
}

std::string System::status_header() const {
  std::stringstream ss;
  for (const Configuration& config : configurations_) {
    ss << config.status_header();
  }
  return ss.str();
}

std::string System::status() const {
  std::stringstream ss;
  for (const Configuration& config : configurations_) {
    ss << config.status();
  }
  return ss.str();
}

void System::synchronize_(const System& system, const Select& perturbed) {
  for (int config = 0; config < num_configurations(); ++config) {
    configurations_[config].synchronize_(system.configuration(config),
      perturbed);
    ASSERT(config == 0, "perturb not implemented for multiple configs");
    // HWH suggest: make perturb a vector, one for each config?
  }
  unoptimized_.synchronize_(system.unoptimized(), perturbed);
  optimized_.synchronize_(system.optimized(), perturbed);
  for (int ref = 0; ref < num_references(); ++ref) {
    references_[ref].synchronize_(system.references()[ref], perturbed);
  }
}

void System::change_volume(const double delta_volume, const argtype& args) {
  Arguments args_(args);
  args_.dont_check();
  const int config = args_.key("configuration").dflt("0").integer();
  const int dimen = args_.key("dimension").dflt("-1").integer();
  configurations_[config].change_volume(delta_volume,
    args_.remove("configuration", args));
  unoptimized_.change_volume(delta_volume, dimen);
  optimized_.change_volume(delta_volume, dimen);
  for (PotentialFactory& ref : references_) {
    ref.change_volume(delta_volume, dimen);
  }
}

}  // namespace feasst
