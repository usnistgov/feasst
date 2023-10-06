#include "utils/include/serialize.h"
#include "configuration/include/neighbor_criteria.h"
#include "monte_carlo/include/trial_select.h"

namespace feasst {

TrialSelect::TrialSelect(argtype args) : TrialSelect(&args) {
  FEASST_CHECK_ALL_USED(args);
}
TrialSelect::TrialSelect(argtype * args) {
  // defaults
  set_ghost(false);

  // parse particle type and group index from args.
  particle_type_ = -1;
  configuration_index_ = integer("configuration_index", args, 0);
  group_index_ = 0;
  if (used("particle_type", *args)) {
    is_particle_type_set_ = true;
    particle_type_ = integer("particle_type", args);
    DEBUG("particle_type " << particle_type_);
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
  DEBUG("is_particle_type_set_ " << is_particle_type_set_);
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
  const int part = anchor_.particle_index(particle_index);
  const int site = anchor_.site_index(particle_index, site_index);
  DEBUG("site " << site);
  return configuration(system).select_particle(part).site(site).position();
}

void TrialSelect::set_ghost(const bool ghost) {
  is_ghost_ = ghost;
  if (is_ghost_) {
    // ASSERT(group_index() == 0, "ghost particles cannot be selected by groups");
    ASSERT(particle_type() >= 0, "ghost particles must be selected by type");
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
  feasst_serialize_version(274, ostr);
  feasst_serialize_fstobj(mobile_original_, ostr);
  feasst_serialize_fstobj(mobile_, ostr);
  feasst_serialize_fstobj(anchor_, ostr);
  feasst_serialize(group_index_, ostr);
  feasst_serialize(particle_type_, ostr);
  feasst_serialize(configuration_index_, ostr);
  feasst_serialize(is_particle_type_set_, ostr);
  feasst_serialize(is_ghost_, ostr);
  feasst_serialize_fstobj(properties_, ostr);
}

TrialSelect::TrialSelect(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 273 && version <= 274, "mismatch version: " << version);
  feasst_deserialize_fstobj(&mobile_original_, istr);
  feasst_deserialize_fstobj(&mobile_, istr);
  feasst_deserialize_fstobj(&anchor_, istr);
  feasst_deserialize(&group_index_, istr);
  feasst_deserialize(&particle_type_, istr);
  if (version >= 274) {
    feasst_deserialize(&configuration_index_, istr);
  }
  feasst_deserialize(&is_particle_type_set_, istr);
  feasst_deserialize(&is_ghost_, istr);
  feasst_deserialize_fstobj(&properties_, istr);
}

void TrialSelect::remove_unphysical_sites(const Configuration& config) {
  //Select unphysical;
  bool resize = false;
  for (int sp_index = 0;
       sp_index < static_cast<int>(mobile_.particle_indices().size());
       ++sp_index) {
    const int p_index = mobile_.particle_indices()[sp_index];
    DEBUG("p_index " << p_index);
    std::vector<int> sites;
    for (const int s_index : mobile_.site_indices(sp_index)) {
      DEBUG("s_index " << s_index);
      if (!config.select_particle(p_index).site(s_index).is_physical()) {
        DEBUG("unphysical");
        sites.push_back(s_index);
      }
    }
    if (sites.size() > 0) {
      resize = true;
      mobile_.remove_sites(p_index, sites);
    }
  }
  if (resize) {
    mobile_.resize_positions();
    mobile_.load_positions(config.particles());
  }
}

void TrialSelect::replace_mobile(const Select& replacement,
    const int sp_index,
    const Configuration& config) {
  bool fast = mobile_.replace_indices(replacement.particle_index(sp_index),
                                      replacement.site_indices(sp_index));
  if (!fast) mobile_.resize_positions();
  mobile_.load_positions(config.particles());
}

bool TrialSelect::select(
    const Select& perturbed,
    System * system,
    Random * random) {
  FATAL("not implemented");
}

const EnergyMap& TrialSelect::map_(const System& system,
    const int neighbor_index) const {
  const int num_neigh = static_cast<int>(system.neighbor_criteria(configuration_index()).size());
  ASSERT(num_neigh > neighbor_index,
    "With " << num_neigh << " NeighborCriteria added to system, the index "
    << neighbor_index << " is out of range.");
  const NeighborCriteria& neighbor_criteria =
    system.neighbor_criteria(neighbor_index, configuration_index());
  DEBUG("ref potential " << neighbor_criteria.reference_potential());
  if (neighbor_criteria.reference_potential() == -1) {
    return system.potentials().potential(neighbor_criteria.potential_index()).visit_model().inner().energy_map();
  }
  const Potential& ref = system.reference(
    neighbor_criteria.reference_potential(),
    neighbor_criteria.potential_index());
  return ref.visit_model().inner().energy_map();
}

void TrialSelect::before_select() {
  mobile_.reset_excluded_and_bond();
  //exclude_energy_ = 0.;
}

void TrialSelect::add_exclude_energy(const double energy) {
  exclude_energy_ = energy;
}

bool TrialSelect::is_isotropic(const System * system) const {
  const Configuration& config = configuration(*system);
  if (config.model_params().index("anisotropic") == -1) {
    return true;
  } else {
    return false;
  }
}

void TrialSelect::set_mobile_original(const System * system) {
  mobile_original_ = mobile_;
  DEBUG("is system isotropic? " << is_isotropic(system));
  if (!is_isotropic(system)) {
    const Configuration& config = configuration(*system);
    for (int select_index = 0;
         select_index < mobile_original_.num_particles();
         ++select_index) {
      const int part_index = mobile_original_.particle_index(select_index);
      for (int select_site = 0;
           select_site < static_cast<int>(mobile_original_.site_indices(select_index).size());
           ++select_site) {
        const int site_index = mobile_original_.site_index(select_index, select_site);
        const Particle& part = config.select_particle(part_index);
        const Site& site = part.site(site_index);
        mobile_original_.set_euler(select_index, select_site, site.euler());
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

}  // namespace feasst
