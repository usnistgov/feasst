#include "utils/include/serialize.h"
#include "math/include/random.h"
#include "chain/include/select_two_sites.h"

namespace feasst {

SelectTwoSites::SelectTwoSites(argtype args) : SelectTwoSites(&args) {
  FEASST_CHECK_ALL_USED(args);
}
SelectTwoSites::SelectTwoSites(argtype * args) : TrialSelect(args) {
  class_name_ = "SelectTwoSites";
  mobile_site_ = integer("mobile_site", args);
  mobile_site2_ = integer("mobile_site2", args);
  ASSERT(mobile_site_ != mobile_site2_, "the mobile site: " << mobile_site_ <<
    " cannot be the same as mobile_site2_: " << mobile_site2_);
}

class MapSelectTwoSites {
 public:
  MapSelectTwoSites() {
    auto obj = MakeSelectTwoSites({{"mobile_site", "1"}, {"mobile_site2", "0"}});
    obj->deserialize_map()["SelectTwoSites"] = obj;
  }
};

static MapSelectTwoSites mapper_ = MapSelectTwoSites();

void SelectTwoSites::precompute(System * system) {
  TrialSelect::precompute(system);
  mobile_.clear();
  mobile_.add_site(0, mobile_site_);
  mobile_.add_site(0, mobile_site2_);
}

bool SelectTwoSites::select(const Select& perturbed,
                             System * system,
                             Random * random) {
  Configuration * config = system->get_configuration();
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
  mobile_.set_particle(0, particle_index);
  mobile_.load_positions(config->particles());
  DEBUG("mobile: " << mobile_.str());
  set_mobile_original(system);
  return true;
}

std::shared_ptr<TrialSelect> SelectTwoSites::create(std::istream& istr) const {
  return std::make_shared<SelectTwoSites>(istr);
}

SelectTwoSites::SelectTwoSites(std::istream& istr) : TrialSelect(istr) {
  // ASSERT(class_name_ == "SelectTwoSites", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 1365, "mismatch version: " << version);
  feasst_deserialize(&mobile_site_, istr);
  feasst_deserialize(&mobile_site2_, istr);
}

void SelectTwoSites::serialize_select_two_sites_(std::ostream& ostr) const {
  serialize_trial_select_(ostr);
  feasst_serialize_version(1365, ostr);
  feasst_serialize(mobile_site_, ostr);
  feasst_serialize(mobile_site2_, ostr);
}

void SelectTwoSites::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_select_two_sites_(ostr);
}

}  // namespace feasst
