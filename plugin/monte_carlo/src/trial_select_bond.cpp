#include "monte_carlo/include/trial_select_bond.h"

namespace feasst {

TrialSelectBond::TrialSelectBond(const argtype& args) : TrialSelect(args) {
  class_name_ = "TrialSelectBond";
  Arguments args_(args);
  args_.dont_check();
  mobile_site_ = args_.key("mobile_site").integer();
  anchor_site_ = args_.key("anchor_site").integer();
}

class MapTrialSelectBond {
 public:
  MapTrialSelectBond() {
    auto obj = MakeTrialSelectBond({{"mobile_site", "1"}, {"anchor_site", "0"}});
    obj->deserialize_map()["TrialSelectBond"] = obj;
  }
};

static MapTrialSelectBond mapper_ = MapTrialSelectBond();

void TrialSelectBond::precompute(System * system) {
  TrialSelect::precompute(system);
  const Particle& part = system->configuration().particle_types().particle(particle_type());
  const int bond_type = part.bond(mobile_site_, anchor_site_).type();
  const Bond& bond = system->configuration().unique_types().particle(particle_type()).bond(bond_type);
  add_property("bond_length", bond.property("length"));
  anchor_.clear();
  anchor_.add_site(0, anchor_site_);
  mobile_.clear();
  mobile_.add_site(0, mobile_site_);
}

bool TrialSelectBond::select(const Select& perturbed,
                             System * system,
                             Random * random) {
  Configuration * config = system->get_configuration();
  int particle_index = -1;
  if (perturbed.num_sites() > 0) {
    particle_index = perturbed.particle_indices().back();
    set_probability(1.);
  } else {
    // select random particle of correct type
    const int group_index = config->particle_type_to_group_create(particle_type());
    const int num = config->num_particles(group_index);
    if (num <= 0) return false;
    const int index = random->uniform(0, num - 1);
    const SelectGroup& select = config->group_select(group_index);
    particle_index = select.particle_index(index);
    set_probability(1./static_cast<double>(num));
  }
  mobile_.set_particle(0, particle_index);
  anchor_.set_particle(0, particle_index);
  mobile_.load_positions(config->particles());
  DEBUG("mobile: " << mobile_.str());
  DEBUG("anchor: " << anchor_.str());
  mobile_original_ = mobile_;
  return true;
}

std::shared_ptr<TrialSelect> TrialSelectBond::create(std::istream& istr) const {
  return std::make_shared<TrialSelectBond>(istr);
}

TrialSelectBond::TrialSelectBond(std::istream& istr)
  : TrialSelect(istr) {
  // ASSERT(class_name_ == "TrialSelectBond", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(235 == version, "mismatch version: " << version);
  feasst_deserialize(&mobile_site_, istr);
  feasst_deserialize(&anchor_site_, istr);
}

void TrialSelectBond::serialize_trial_select_bond_(std::ostream& ostr) const {
  serialize_trial_select_(ostr);
  feasst_serialize_version(235, ostr);
  feasst_serialize(mobile_site_, ostr);
  feasst_serialize(anchor_site_, ostr);
}

void TrialSelectBond::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_select_bond_(ostr);
}

}  // namespace feasst
