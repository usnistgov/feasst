#include <cmath>  // isinf and isnan
#include "utils/include/io.h"
#include "utils/include/arguments.h"
#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "configuration/include/domain.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/select.h"
#include "configuration/include/neighbor_criteria.h"
#include "configuration/include/configuration.h"
#include "system/include/bond_visitor.h"
#include "system/include/thermo_params.h"
#include "system/include/system.h"

namespace feasst {

void System::add(std::shared_ptr<Configuration> configuration) {
  configurations_.push_back(configuration);
  bonds_.push_back(std::make_shared<BondVisitor>());
  unoptimized_.push_back(PotentialFactory());
  optimized_.push_back(PotentialFactory());
}

const Configuration& System::configuration(const int config) const {
  return *configurations_[config];
}

Configuration* System::get_configuration(const int config) {
  return configurations_[config].get();
}

int System::dimension(const int config) const {
  const int dim = configurations_[0]->dimension();
  for (auto config : configurations_) {
    ASSERT(dim == config->dimension(), "dimensions of configs do not match");
  }
  return dim;
}

void System::add_to_unoptimized(std::shared_ptr<Potential> potential,
    const int config) {
  unoptimized_[config].add(potential);
  unoptimized_[config].precompute(unoptimized_[config].num() - 1,
                                               configurations_[config].get());
}

void System::set_unoptimized(const int index,
    std::shared_ptr<Potential> potential,
    const int config) {
  unoptimized_[config].set(index, potential);
  unoptimized_[config].precompute(index, configurations_[config].get());
}

void System::add_to_optimized(std::shared_ptr<Potential> potential,
    const int config) {
  is_optimized_ = true;
  optimized_[config].add(potential);
  optimized_[config].precompute(optimized_[config].num() - 1,
                                configurations_[config].get());
}

PotentialFactory * System::reference_(const int index, const int config) {
  ASSERT(index < static_cast<int>(references_[config].size()),
    "An unrecognized reference potential: " << index << " was requested. "
    << "But there are only " << references_[config].size() << " reference potentials. "
    << "Perhaps check if the desired \"RefPotential\" is in the input script.");
  return &references_[config][index];
}

void System::add_to_reference(std::shared_ptr<Potential> ref, const int index,
    const int config) {
  if (config == static_cast<int>(references_.size())) {
    references_.push_back(std::vector<PotentialFactory>());
  } else if (index > static_cast<int>(references_.size())) {
    FATAL("Add references in order of the configuration.");
  }
  if (index == static_cast<int>(references_[config].size())) {
    references_[config].push_back(PotentialFactory());
  } else if (index > static_cast<int>(references_[config].size())) {
    FATAL("references must be added in order");
  }
  reference_(index, config)->add(ref);
  reference_(index, config)->precompute(reference_(index, config)->num() - 1,
                                        configurations_[config].get());
}

const Potential& System::reference(const int ref,
    const int potential,
    const int config) const {
  ASSERT(ref < static_cast<int>(references_[config].size()), "reference potential: "
    << ref << " >= number of references: " << references_[config].size()
    << ". Double check the reference_index.");
  return references_[config][ref].potential(potential);
}

void System::precompute() {
  for (int config_index = 0; config_index < num_configurations(); ++config_index) {
    Configuration * config = configurations_[config_index].get();
    unoptimized_[config_index].precompute(config);
    if (is_optimized_) {
      optimized_[config_index].precompute(config);
    }
    if (num_references(config_index) > 0) {
      for (PotentialFactory& ref : references_[config_index]) {
        ref.precompute(config);
      }
    }
  }
}

double System::unoptimized_energy(const int config) {
  ASSERT(config < num_configurations(),
    "config: " << config << " >= num_configurations: " << num_configurations());
  const double en = unoptimized_[config].energy(configurations_[config].get());
  ASSERT(!std::isinf(en) && !std::isnan(en),
    "Energy(" << en << ") is infinite or not "
    << "a number. Are particles on top of each other?");
  unoptimized_[config].finalize(configurations_[config]->selection_of_all(),
    configurations_[config].get());
  ref_used_last_ = -1;
  DEBUG("ref_used_last_ " << ref_used_last_);
  bonds_[config]->compute_all(*configurations_[config]);
  return en + bonds_[config]->energy();
}

PotentialFactory * System::potentials_(const int config) {
  if (is_optimized_) {
    return &optimized_[config];
  }
  return &unoptimized_[config];
}

double System::energy(const int config) {
  const double en = potentials_(config)->energy(configurations_[config].get());
  finalize(config);
  ref_used_last_ = -1;
  DEBUG("ref_used_last_ " << ref_used_last_);
  bonds_[config]->compute_all(*configurations_[config]);
  DEBUG("bond en " << bonds_[config]->energy());
  return en + bonds_[config]->energy();
}

double System::perturbed_energy(const Select& select, const int config) {
  ref_used_last_ = -1;
  DEBUG("ref_used_last_ " << ref_used_last_);
  double en = potentials_(config)->select_energy(select, configurations_[config].get());
  bonds_[config]->compute_all(select, *configurations_[config]);
  const double bond_en = bonds_[config]->energy();
  DEBUG("bond en " << bonds_[config]->energy());
  ASSERT(!std::isinf(en), "en: " << en << " is inf.");
  ASSERT(!std::isnan(en), "en: " << en << " is nan.");
  ASSERT(!std::isinf(bond_en), "bond_en: " << bond_en << " is inf.");
  ASSERT(!std::isnan(bond_en), "bond_en: " << bond_en << " is nan.");
  return en + bond_en;
}

double System::reference_energy(const int ref, const int config) {
  ref_used_last_ = ref;
  DEBUG("ref_used_last_ " << ref_used_last_);
  return reference_(ref, config)->energy(configurations_[config].get());
}

double System::reference_energy(const Select& select,
    const int ref,
    const int config) {
  ref_used_last_ = ref;
  DEBUG("ref_used_last_ " << ref_used_last_);
  ASSERT(ref < num_references(), "Asked for reference: " << ref <<
    ", but there are only " << num_references() << " RefPotentials.");
  return reference_(ref, config)->select_energy(select,
                                                configurations_[config].get());
}

void System::serialize(std::ostream& sstr) const {
  feasst_serialize_version(7349, sstr);
  feasst_serialize(configurations_, sstr);
  feasst_serialize(bonds_, sstr);
  feasst_serialize_fstobj(unoptimized_, sstr);
  feasst_serialize_fstobj(optimized_, sstr);
  feasst_serialize(is_optimized_, sstr);
  feasst_serialize_fstobj(references_, sstr);
  feasst_serialize(thermo_params_, sstr);
  feasst_serialize_endcap("System", sstr);
  DEBUG("size: " << sstr.tellp());
}

System::System(std::istream& sstr) {
  const int version = feasst_deserialize_version(sstr);
  ASSERT(version == 7349, "unrecognized verison: " << version);
//  HWH for unknown reasons, this function template does not work.
  //feasst_deserialize(&configurations_, sstr);
  {
    int dim1;
    sstr >> dim1;
    configurations_.resize(dim1);
    for (int index = 0; index < dim1; ++index) {
      //feasst_deserialize((*vector)[index], istr);
      int existing;
      sstr >> existing;
      if (existing != 0) {
        configurations_[index] = std::make_shared<Configuration>(sstr);
      }
    }
  }
//  HWH for unknown reasons, this function template does not work.
  //feasst_deserialize(&bonds_, sstr);
  {
    int dim1;
    sstr >> dim1;
    bonds_.resize(dim1);
    for (int index = 0; index < dim1; ++index) {
      //feasst_deserialize((*vector)[index], istr);
      int existing;
      sstr >> existing;
      if (existing != 0) {
        bonds_[index] = std::make_shared<BondVisitor>(sstr);
      }
    }
  }
  feasst_deserialize_fstobj(&unoptimized_, sstr);
  feasst_deserialize_fstobj(&optimized_, sstr);
  feasst_deserialize(&is_optimized_, sstr);
  feasst_deserialize_fstobj(&references_, sstr);
  // HWH for unknown reasons, this function template does not work.
  //feasst_deserialize(thermo_params_, sstr);
  { int existing;
    sstr >> existing;
    if (existing != 0) {
      thermo_params_ = std::make_shared<ThermoParams>(sstr);
    }
  }
  feasst_deserialize_endcap("System", sstr);
}

const PotentialFactory& System::potentials(const int config) const {
  if (is_optimized_) {
    return optimized_[config];
  }
  return unoptimized_[config];
}

void System::load_cache(const bool load) {
  for (int config = 0; config < num_configurations(); ++config) {
    unoptimized_[config].load_cache(load);
    optimized_[config].load_cache(load);
    if (num_references(config) > 0) {
      for (PotentialFactory& ref : references_[config]) {
        ref.load_cache(load);
      }
    }
  }
}

void System::unload_cache(const System& system) {
  for (int config = 0; config < num_configurations(); ++config) {
    unoptimized_[config].unload_cache(system.unoptimized_[config]);
    optimized_[config].unload_cache(system.optimized_[config]);
    if (num_references(config) > 0) {
      ASSERT(references_[config].size() == system.references_[config].size(),
             "size mismatch");
      for (int iref = 0; iref < static_cast<int>(references_[config].size()); ++iref) {
        references_[config][iref].unload_cache(system.references_[config][iref]);
      }
    }
  }
}

void System::finalize(const Select& select, const int config) {
  DEBUG("finalizing " << select.str() << " in config:" << config);
  if (select.trial_state() == 2) {
    // finalize removal
    configurations_[config]->remove_particles(select);
  }
  unoptimized_[config].finalize(select, configurations_[config].get());
  optimized_[config].finalize(select, configurations_[config].get());
  if (num_references(config) > 0) {
    for (PotentialFactory& ref : references_[config]) {
      ref.finalize(select, configurations_[config].get());
    }
  }
  for (int iconf = 0; iconf < num_configurations(); ++iconf) {
    DEBUG("number particles in conf " << iconf << ": " << configuration(iconf).num_particles());
  }
}

void System::revert(const Select& select, const int config) {
  DEBUG("ref_used_last_ " << ref_used_last_);
  if (ref_used_last_ != -1) {
    references_[config][ref_used_last_].revert(select);
  } else {
    potentials_(config)->revert(select);
  }
}

void System::check(const int config) const {
  unoptimized_[config].check(*configurations_[config]);
  optimized_[config].check(*configurations_[config]);
  if (num_references(config) > 0) {
    for (const PotentialFactory& ref : references_[config]) {
      ref.check(*configurations_[config]);
    }
  }
}

std::string System::status_header() const {
  std::stringstream ss;
  if (num_configurations() == 1) {
    ss << configuration(0).status_header();
  } else {
    for (int iconf = 0; iconf < num_configurations(); ++iconf) {
      const Configuration& config = *configurations_[iconf];
      ss << config.status_header("_config" + str(iconf));
    }
  }
  ss << ",beta";
  return ss.str();
}

std::string System::status() const {
  std::stringstream ss;
  for (auto config : configurations_) {
    ss << config->status();
  }
  ss << "," << thermo_params().beta();
  return ss.str();
}

  // HWH suggest: make perturb a vector, one for each config?
void System::synchronize_(const System& system, const std::vector<std::shared_ptr<Select> >& perturbed) {
  for (int config = 0; config < num_configurations(); ++config) {
    configurations_[config]->synchronize_(system.configuration(config),
      *perturbed[config]);
    unoptimized_[config].synchronize_(system.unoptimized(), *perturbed[config]);
    optimized_[config].synchronize_(system.optimized(), *perturbed[config]);
    for (int ref = 0; ref < num_references(config); ++ref) {
      references_[config][ref].synchronize_(system.references()[config][ref], *perturbed[config]);
    }
  }
}

void System::change_volume(const double delta_volume, argtype * args) {
  const int config = integer("configuration", args, 0);
  const int dimen = integer("dimension", args, -1);
  args->insert({"dimension", str(dimen)});
  configurations_[config]->change_volume(delta_volume, args);
  unoptimized_[config].change_volume(delta_volume, dimen, configurations_[config].get());
  optimized_[config].change_volume(delta_volume, dimen, configurations_[config].get());
  if (num_references(config) > 0) {
    for (PotentialFactory& ref : references_[config]) {
      ref.change_volume(delta_volume, dimen, configurations_[config].get());
    }
  }
  delta_volume_previous_ = delta_volume;
  // HWH check cutoff is less than half box
  // HWH or don't change volume at all?
}
void System::change_volume(const double delta_volume, argtype args) {
  change_volume(delta_volume, &args);
  feasst_check_all_used(args);
}

//double System::constrained_volume_change(argtype * args) {
//  change_volume(-delta_volume_previous_, args);
//  return -delta_volume_previous_;
//}
//double System::constrained_volume_change(argtype args) {
//  const double delta_volume = constrained_volume_change(&args);
//  feasst_check_all_used(args);
//  return delta_volume;
//}

int System::num_references(const int config) const {
  if (static_cast<int>(references_.size()) > 0) {
    return static_cast<int>(references_[config].size());
  }
  return 0;
}

void System::set(std::shared_ptr<ThermoParams> thermo_params) {
  thermo_params_ = thermo_params;
}

const ThermoParams& System::thermo_params() const {
  ASSERT(thermo_params_, "must set ThermoParams first.");
  return const_cast<ThermoParams&>(*thermo_params_);
}

void System::remove_opt_overlap() {
  for (int config = 0; config < num_configurations(); ++config) {
    unoptimized_[config].remove_opt_overlap();
    optimized_[config].remove_opt_overlap();
    if (num_references(config) > 0) {
      for (PotentialFactory& ref : references_[config]) {
        ref.remove_opt_overlap();
      }
    }
  }
}

const NeighborCriteria& System::neighbor_criteria(const int index,
    const int config) const {
  return configurations_[config]->neighbor_criteria(index);
}

const std::vector<std::shared_ptr<NeighborCriteria> >& System::neighbor_criteria(
    const int config) const {
  return configurations_[config]->neighbor_criteria();
}

NeighborCriteria * System::get_neighbor_criteria(const int index,
                                                 const int config) {
  return configurations_[config]->get_neighbor_criteria(index);
}

double System::total_volume() const {
  double volume = 0.;
  for (auto config : configurations_) {
    volume += config->domain().volume();
  }
  return volume;
}

const ThermoParams * System::thermo_params_ptr_() const {
  return thermo_params_.get();
}

void System::set_beta(const double beta) { thermo_params_->set_beta(beta); }

int System::num_configurations() const {
  return static_cast<int>(configurations_.size());
}

void System::add(std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const int config) {
  configurations_[config]->add(neighbor_criteria);
}

void System::finalize(const int config) {
  finalize(configurations_[config]->selection_of_all(), config);
}

double System::initialize(const int config) {
  precompute();
  const double en = unoptimized_energy(config);
  energy(config);
  for (int ref = 0; ref < num_references(config); ++ref) {
    reference_energy(ref, config);
  }
  return en;
}

}  // namespace feasst
