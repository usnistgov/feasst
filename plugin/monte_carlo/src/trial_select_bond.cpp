#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "math/include/random.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/configuration.h"
#include "monte_carlo/include/trial_select_bond.h"

namespace feasst {

TrialSelectBond::TrialSelectBond(argtype args) : TrialSelectBond(&args) {
  feasst_check_all_used(args);
}
TrialSelectBond::TrialSelectBond(argtype * args) : TrialSelect(args) {
  class_name_ = "TrialSelectBond";
  mobile_site_name_ = str("mobile_site", args);
  anchor_site_name_ = str("anchor_site", args);
  ASSERT(mobile_site_name_ != anchor_site_name_, "the mobile site: " <<
    mobile_site_name_ << " cannot be the same as the anchor site: " <<
    anchor_site_name_);
  ignore_bond_ = boolean("ignore_bond", args, false);
}

FEASST_MAPPER(TrialSelectBond,
  argtype({{"mobile_site", "1"}, {"anchor_site", "0"}}));

void TrialSelectBond::precompute(System * system) {
  TrialSelect::precompute(system);
  const Configuration& conf = configuration(*system);
  const int mobile_site = conf.site_name_to_index(mobile_site_name_);
  const int anchor_site = conf.site_name_to_index(anchor_site_name_);
  if (!ignore_bond_) {
    const Particle& part = conf.particle_types().particle(particle_type());
    const int bond_type = part.bond(mobile_site, anchor_site).type();
    add_or_set_property("bond_type", bond_type);
  }
  anchor_->clear();
  anchor_->add_site(0, anchor_site);
  get_mobile()->clear();
  get_mobile()->add_site(0, mobile_site);
}

bool TrialSelectBond::select(const Select& perturbed,
                             System * system,
                             Random * random,
                             TrialSelect * previous_select) {
  Configuration * config = get_configuration(system);
  int particle_index = -1;
  if (perturbed.num_sites() > 0) {
    particle_index = perturbed.particle_indices().back();
    set_probability_(1.);
  } else {
    // select random particle of correct type
    const int group_index = config->particle_type_to_group_create(particle_type());
    const int num = config->num_particles(group_index);
    if (num <= 0) return false;
    const int index = random->uniform(0, num - 1);
    const Select& select = config->group_select(group_index);
    particle_index = select.particle_index(index);
    set_probability_(1./static_cast<double>(num));
  }
  get_mobile()->set_particle(0, particle_index);
  get_anchor()->set_particle(0, particle_index);
  get_mobile()->load_positions(config->particles());
  DEBUG("mobile: " << mobile().str());
  DEBUG("anchor: " << anchor().str());
  set_mobile_original(system);
  return true;
}

std::shared_ptr<TrialSelect> TrialSelectBond::create(std::istream& istr) const {
  return std::make_shared<TrialSelectBond>(istr);
}

TrialSelectBond::TrialSelectBond(std::istream& istr)
  : TrialSelect(istr) {
  // ASSERT(class_name_ == "TrialSelectBond", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 235 && version <= 237, "mismatch version: " << version);
  if (version <= 236) {
    WARN("Restart versions may be incompatible");
    int anchor_site;
    feasst_deserialize(&anchor_site, istr);
  }
  if (version >= 237) {
    feasst_deserialize(&mobile_site_name_, istr);
    feasst_deserialize(&anchor_site_name_, istr);
  }
  if (version >= 236) {
    feasst_deserialize(&ignore_bond_, istr);
  }
}

void TrialSelectBond::serialize_trial_select_bond_(std::ostream& ostr) const {
  serialize_trial_select_(ostr);
  feasst_serialize_version(237, ostr);
  feasst_serialize(mobile_site_name_, ostr);
  feasst_serialize(anchor_site_name_, ostr);
  feasst_serialize(ignore_bond_, ostr);
}

void TrialSelectBond::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_select_bond_(ostr);
}

}  // namespace feasst
