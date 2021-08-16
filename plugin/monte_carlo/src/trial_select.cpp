#include "utils/include/serialize.h"
#include "system/include/neighbor_criteria.h"
#include "monte_carlo/include/trial_select.h"

namespace feasst {

TrialSelect::TrialSelect(argtype args) : TrialSelect(&args) {
  check_all_used(args);
}
TrialSelect::TrialSelect(argtype * args) {
  // defaults
  set_ghost(false);

  // parse particle type and group index from args.
  particle_type_ = -1;
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
    group_index_ = system->get_configuration()->particle_type_to_group_create(
      particle_type_);
    DEBUG("group_index_ " << group_index_);
  }
}

const Position& TrialSelect::anchor_position(const int particle_index,
    const int site_index,
    const System& system) {
  const int part = anchor_.particle_index(particle_index);
  const int site = anchor_.site_index(particle_index, site_index);
  DEBUG("site " << site);
  return system.configuration().select_particle(part).site(site).position();
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
  feasst_serialize_version(273, ostr);
  feasst_serialize_fstobj(mobile_original_, ostr);
  feasst_serialize_fstobj(mobile_, ostr);
  feasst_serialize_fstobj(anchor_, ostr);
  feasst_serialize(group_index_, ostr);
  feasst_serialize(particle_type_, ostr);
  feasst_serialize(is_particle_type_set_, ostr);
  feasst_serialize(is_ghost_, ostr);
  feasst_serialize_fstobj(properties_, ostr);
}

TrialSelect::TrialSelect(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(273 == version, "mismatch version: " << version);
  feasst_deserialize_fstobj(&mobile_original_, istr);
  feasst_deserialize_fstobj(&mobile_, istr);
  feasst_deserialize_fstobj(&anchor_, istr);
  feasst_deserialize(&group_index_, istr);
  feasst_deserialize(&particle_type_, istr);
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
  DEBUG("here");
  const int num_neigh = static_cast<int>(system.neighbor_criteria().size());
  ASSERT(num_neigh > neighbor_index,
    "With " << num_neigh << " NeighborCriteria added to system, the index "
    << neighbor_index << " is out of range.");
  const NeighborCriteria& neighbor_criteria =
    system.neighbor_criteria(neighbor_index);
  DEBUG("here");
  DEBUG("ref potential " << neighbor_criteria.reference_potential());
  if (neighbor_criteria.reference_potential() == -1) {
    return system.potentials().potential(neighbor_criteria.potential_index()).visit_model().inner().energy_map();
  }
  return system.reference(neighbor_criteria.reference_potential(),
                          neighbor_criteria.potential_index()
                         ).visit_model().inner().energy_map();
}

}  // namespace feasst
