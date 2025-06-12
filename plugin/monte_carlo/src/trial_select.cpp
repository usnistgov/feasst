#include "utils/include/serialize_extra.h"
#include "utils/include/arguments.h"
#include "math/include/accumulator.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/configuration.h"
#include "configuration/include/neighbor_criteria.h"
#include "configuration/include/model_params.h"
#include "configuration/include/properties.h"
#include "configuration/include/select.h"
#include "system/include/potential.h"
#include "system/include/visit_model.h"
#include "system/include/visit_model_inner.h"
#include "system/include/system.h"
#include "monte_carlo/include/trial_select.h"

namespace feasst {

TrialSelect::TrialSelect(argtype args) : TrialSelect(&args) {
  feasst_check_all_used(args);
}
TrialSelect::TrialSelect(argtype * args) {
  properties_ = std::make_shared<Properties>();
  mobile_original_ = std::make_shared<Select>();
  mobile_ = std::make_shared<Select>();
  anchor_ = std::make_shared<Select>();
  // defaults
  set_ghost(false);

  // parse particle type and group index from args.
  particle_type_ = -1;
//  if (used("configuration_index", *args)) {
//    WARN("Deprecated TrialSelect::configuration_index->config (see Configuration::name)");
//  }
  configuration_index_ = integer("configuration_index", args, 0);
  config_ = str("config", args, "");
  group_index_ = 0;
  if (used("particle_type", *args)) {
    particle_type_name_ = str("particle_type", args);
    DEBUG("particle_type " << particle_type_name_);
    ASSERT(!used("group_index", *args),
      "cant specify both particle type and group index");
  } else {
    if (used("group_index", *args)) {
      group_index_ = integer("group_index", args);
      ASSERT(!used("group", *args),
        "cant specify both group_index and group name");
    } else {
      if (used("group", *args)) {
        group_ = str("group", args);
      }
    }
  }

  set_probability_();
}

int TrialSelect::particle_type() const {
  ASSERT(is_particle_type_set_, "particle type not specified");
  return particle_type_;
}

void TrialSelect::precompute(System * system) {
  if (!config_.empty()) {
    configuration_index_ = system->configuration_index(config_);
  }
  aniso_index_ = system->configuration().model_params().index("anisotropic");
  DEBUG("is_particle_type_set_ " << is_particle_type_set_);
  if (!particle_type_name_.empty()) {
    is_particle_type_set_ = true;
    const Configuration& config = configuration(*system);
    particle_type_ = config.particle_name_to_type(particle_type_name_);
  }
  if (is_particle_type_set_) {
    DEBUG("particle_type " << particle_type_);
    group_index_ = get_configuration(system)->particle_type_to_group_create(
      particle_type_);
    DEBUG("group_index_ " << group_index_);
  } else if (!group_.empty()) {
    group_index_ = configuration(*system).group_index(group_);
  }
}

const Position& TrialSelect::anchor_position(const int particle_index,
    const int site_index,
    const System& system) const {
  const int part = anchor_->particle_index(particle_index);
  const int site = anchor_->site_index(particle_index, site_index);
  DEBUG("site " << site);
  return configuration(system).select_particle(part).site(site).position();
}

void TrialSelect::set_ghost(const bool ghost) {
  is_ghost_ = ghost;
  if (is_ghost_) {
    // ASSERT(group_index() == 0, "ghost particles cannot be selected by groups");
    ASSERT(!particle_type_name().empty(), "ghost particles must be selected by type");
  }
}

std::map<std::string, std::shared_ptr<TrialSelect> >& TrialSelect::deserialize_map() {
  static std::map<std::string, std::shared_ptr<TrialSelect> >* ans =
     new std::map<std::string, std::shared_ptr<TrialSelect> >();
  return *ans;
}

void TrialSelect::serialize(std::ostream& ostr) const { FATAL("not implemented"); }

std::shared_ptr<TrialSelect> TrialSelect::create(std::istream& istr) const {
  FATAL("not implemented");
}

std::shared_ptr<TrialSelect> TrialSelect::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

void TrialSelect::serialize_trial_select_(std::ostream& ostr) const {
  feasst_serialize_version(277, ostr);
  feasst_serialize(mobile_original_, ostr);
  feasst_serialize(mobile_, ostr);
  feasst_serialize(anchor_, ostr);
  feasst_serialize(printable_, ostr);
  feasst_serialize(group_index_, ostr);
  feasst_serialize(particle_type_, ostr);
  feasst_serialize(particle_type_name_, ostr);
  feasst_serialize(configuration_index_, ostr);
  feasst_serialize(config_, ostr);
  feasst_serialize(is_particle_type_set_, ostr);
  feasst_serialize(is_ghost_, ostr);
  feasst_serialize(properties_, ostr);
  feasst_serialize(aniso_index_, ostr);
}

TrialSelect::TrialSelect(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 273 && version <= 277, "mismatch version: " << version);
//  feasst_deserialize(mobile_original_, istr);
// HWH for unknown reasons, this function template does not work.
  {
    int existing;
    istr >> existing;
    if (existing != 0) {
      mobile_original_ = std::make_shared<Select>(istr);
    }
  }
//  feasst_deserialize(mobile_, istr);
// HWH for unknown reasons, this function template does not work.
  {
    int existing;
    istr >> existing;
    if (existing != 0) {
      mobile_ = std::make_shared<Select>(istr);
    }
  }
//  feasst_deserialize(anchor_, istr);
// HWH for unknown reasons, this function template does not work.
  {
    int existing;
    istr >> existing;
    if (existing != 0) {
      anchor_ = std::make_shared<Select>(istr);
    }
  }
  if (version >= 275) {
    feasst_deserialize(&printable_, istr);
  }
  feasst_deserialize(&group_index_, istr);
  feasst_deserialize(&particle_type_, istr);
  if (version >= 276) {
    feasst_deserialize(&particle_type_name_, istr);
  }
  if (version >= 274) {
    feasst_deserialize(&configuration_index_, istr);
  }
  if (version >= 277) {
    feasst_deserialize(&config_, istr);
  }
  feasst_deserialize(&is_particle_type_set_, istr);
  feasst_deserialize(&is_ghost_, istr);
//  feasst_deserialize(properties_, istr);
// HWH for unknown reasons, this function template does not work.
  {
    int existing;
    istr >> existing;
    if (existing != 0) {
       properties_ = std::make_shared<Properties>(istr);
    }
  }
  feasst_deserialize(&aniso_index_, istr);
}

void TrialSelect::remove_unphysical_sites(const Configuration& config) {
  //Select unphysical;
  bool resize = false;
  for (int sp_index = 0;
       sp_index < static_cast<int>(mobile_->particle_indices().size());
       ++sp_index) {
    const int p_index = mobile_->particle_indices()[sp_index];
    DEBUG("p_index " << p_index);
    std::vector<int> sites;
    for (const int s_index : mobile_->site_indices(sp_index)) {
      DEBUG("s_index " << s_index);
      if (!config.select_particle(p_index).site(s_index).is_physical()) {
        DEBUG("unphysical");
        sites.push_back(s_index);
      }
    }
    if (sites.size() > 0) {
      resize = true;
      mobile_->remove_sites(p_index, sites);
    }
  }
  if (resize) {
    mobile_->resize_positions();
    mobile_->load_positions(config.particles());
  }
}

void TrialSelect::replace_mobile(const Select& replacement,
    const int sp_index,
    const Configuration& config) {
  bool fast = mobile_->replace_indices(replacement.particle_index(sp_index),
                                       replacement.site_indices(sp_index));
  if (!fast) mobile_->resize_positions();
  mobile_->load_positions(config.particles());
}

bool TrialSelect::select(
    const Select& perturbed,
    System * system,
    Random * random,
    TrialSelect * previous_select) {
  FATAL("not implemented");
}

const EnergyMap& TrialSelect::map_(const System& system,
    const int neighbor_index) const {
  const int iconf = configuration_index();
  const int num_neigh = static_cast<int>(system.neighbor_criteria(iconf).size());
  ASSERT(num_neigh > neighbor_index,
    "With " << num_neigh << " NeighborCriteria added to system, the index "
    << neighbor_index << " is out of range.");
  const NeighborCriteria& neighbor_criteria =
    system.neighbor_criteria(neighbor_index, iconf);
  int iref = neighbor_criteria.reference_potential();
  // HWH optimize by adding to NeighborCriteria.name_to_index given config names or other precompute
  if (!neighbor_criteria.ref().empty()) {
    iref = system.reference_index(iconf, neighbor_criteria.ref());
  }
  const int ipot = neighbor_criteria.potential_index();
  DEBUG("iref:" << iref);
  if (iref == -1) {
    return system.potentials().potential(ipot).visit_model().inner().energy_map();
  }
  const Potential& ref = system.reference(iref, ipot);
  return ref.visit_model().inner().energy_map();
}

void TrialSelect::before_select() {
  mobile_->reset_excluded_and_bond();
  //exclude_energy_ = 0.;
}

void TrialSelect::add_exclude_energy(const double energy) {
  exclude_energy_ = energy;
}

bool TrialSelect::is_isotropic(const System * system) const {
  if (aniso_index_ == -1) {
    return true;
  } else {
    return false;
  }
}

void TrialSelect::set_mobile_original(const System * system) {
  mobile_original_ = std::make_shared<Select>(*mobile_);
  DEBUG("is system isotropic? " << is_isotropic(system));
  if (!is_isotropic(system)) {
    const Configuration& config = configuration(*system);
    for (int select_index = 0;
         select_index < mobile_original_->num_particles();
         ++select_index) {
      const int part_index = mobile_original_->particle_index(select_index);
      for (int select_site = 0;
           select_site < static_cast<int>(mobile_original_->site_indices(select_index).size());
           ++select_site) {
        const int site_index = mobile_original_->site_index(select_index, select_site);
        const Particle& part = config.select_particle(part_index);
        const Site& site = part.site(site_index);
        mobile_original_->set_euler(select_index, select_site, site.euler());
        DEBUG("original Euler(" << part_index << "," << site_index << ") " << site.euler().str());
      }
    }
  }
}

bool TrialSelect::are_constraints_satisfied(const int old,
    const System& system) const {
  return true;
}

void TrialSelect::set_configuration_index(const int config) {
  configuration_index_ = config;
}

void TrialSelect::set_config(const std::string& config) {
  config_ = config;
}

const Configuration& TrialSelect::configuration(const System& system) const {
  return system.configuration(configuration_index_);
}

Configuration * TrialSelect::get_configuration(System * system) const {
  return system->get_configuration(configuration_index_);
}

bool TrialSelect::sel(System * system, Random * random) {
  if (!empty_) {
    empty_ = std::make_shared<Select>();
  }
  return select(*empty_, system, random, NULL);
}

const std::map<std::string, std::shared_ptr<Accumulator> >& TrialSelect::printable() const {
  return printable_;
}

const Accumulator& TrialSelect::printable(const std::string str) const {
  return const_cast<const Accumulator&>(*printable_.at(str));
}

double TrialSelect::property(const std::string name) const {
  return properties_->value(name);
}

bool TrialSelect::has_property(const std::string name) const {
  return properties_->has(name);
}

void TrialSelect::add_or_set_property(const std::string name, const double value) {
  properties_->add_or_set(name, value);
}

const Select& TrialSelect::anchor() const { return *anchor_; }

Select * TrialSelect::get_anchor() { return anchor_.get(); }

const Select& TrialSelect::mobile() const { return *mobile_; }

Select * TrialSelect::get_mobile() { return mobile_.get(); }

void TrialSelect::set_mobile(const Select& mobile) { mobile_ = std::make_shared<Select>(mobile); }

const Select& TrialSelect::mobile_original() const { return *mobile_original_; }

void TrialSelect::set_trial_state(const int state) { mobile_->set_trial_state(state); }

void TrialSelect::reset_mobile() { mobile_ = std::make_shared<Select>(*mobile_original_); }

}  // namespace feasst
