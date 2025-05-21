#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "math/include/utils_math.h"
#include "math/include/random.h"
#include "configuration/include/select.h"
#include "configuration/include/configuration.h"
#include "chain/include/select_crankshaft_small.h"

namespace feasst {

FEASST_MAPPER(SelectCrankshaftSmall,
  argtype({{"site", "0"}, {"anchor_site0", "1"}, {"anchor_site1", "2"}}));

std::shared_ptr<TrialSelect> SelectCrankshaftSmall::create(std::istream& istr) const {
  return std::make_shared<SelectCrankshaftSmall>(istr);
}

SelectCrankshaftSmall::SelectCrankshaftSmall(std::istream& istr)
  : TrialSelectParticle(istr) {
  // ASSERT(class_name_ == "SelectCrankshaftSmall", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(version >= 9287 && version <= 9288, "mismatch version: " << version);
  if (version >= 9288) {
    feasst_deserialize(&site_names_, istr);
  }
}

void SelectCrankshaftSmall::serialize_select_crankshaft_small_(std::ostream& ostr) const {
  serialize_trial_select_particle_(ostr);
  feasst_serialize_version(9288, ostr);
  feasst_serialize(site_names_, ostr);
}

void SelectCrankshaftSmall::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_select_crankshaft_small_(ostr);
}

SelectCrankshaftSmall::SelectCrankshaftSmall(argtype args) : SelectCrankshaftSmall(&args) {
  feasst_check_all_used(args);
}
SelectCrankshaftSmall::SelectCrankshaftSmall(argtype * args) : TrialSelectParticle(args) {
  class_name_ = "SelectCrankshaftSmall";
  DEBUG("parse mobile sites");
  std::string start = "site";
  std::stringstream key;
  int num = 1;
  key << start << num;
  DEBUG("key " << key.str());
  DEBUG("args " << str(*args));
  while (used(key.str(), *args)) {
    site_names_.push_back(str(key.str(), args));
    ++num;
    ASSERT(num < 1e8, "num(" << num << ") is very high. Infinite loop?");
    key.str("");
    key << start << num;
  }
  get_anchor()->clear();
  get_anchor()->add_site(0, integer("anchor_site0", args));
  get_anchor()->add_site(0, integer("anchor_site1", args));
  DEBUG("site names " << feasst_str(site_names_));
  DEBUG("args " << str(*args));
}

void SelectCrankshaftSmall::precompute(System * system) {
  TrialSelectParticle::precompute(system);
  const Configuration& conf = configuration(*system);
  DEBUG("site names " << feasst_str(site_names_));
  for (const std::string& name : site_names_) {
    get_mobile()->add_site(0, conf.site_name_to_index(name));
  }
  ASSERT(mobile().num_sites() > 0, "mobile().num_sites() == " <<
    mobile().num_sites() << ". site argument required");
}

bool SelectCrankshaftSmall::select(const Select& perturbed,
    System * system,
    Random * random,
    TrialSelect * previous_select) {
  Configuration * config = get_configuration(system);
  const int group_index = config->particle_type_to_group_create(particle_type());
  const int num = config->num_particles(group_index);
  if (num <= 0) return false;
  const int index = random->uniform(0, num - 1);
  const Select& select = config->group_select(group_index);
  const int particle_index = select.particle_index(index);
  set_probability_(1./static_cast<double>(num));
  get_mobile()->set_particle(0, particle_index);
  get_anchor()->set_particle(0, particle_index);
  get_mobile()->load_positions(config->particles());
  set_mobile_original(system);
  return true;
}

}  // namespace feasst
