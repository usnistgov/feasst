#include "utils/include/serialize.h"
#include "utils/include/utils.h"  // find_in_list
#include "math/include/random.h"
#include "monte_carlo/include/trial_select_all.h"

namespace feasst {

TrialSelectAll::TrialSelectAll(argtype * args) : TrialSelect(args) {
  class_name_ = "TrialSelectAll";
}
TrialSelectAll::TrialSelectAll(argtype args) : TrialSelectAll(&args) {
  FEASST_CHECK_ALL_USED(args);
}

class MapTrialSelectAll {
 public:
  MapTrialSelectAll() {
    auto obj = MakeTrialSelectAll();
    obj->deserialize_map()["TrialSelectAll"] = obj;
  }
};

static MapTrialSelectAll mapper_ = MapTrialSelectAll();

bool TrialSelectAll::select(const Select& perturbed,
                                 System* system,
                                 Random * random) {
  set_probability_(1.);
  const Configuration& config = system->configuration();
  set_mobile(config.selection_of_all());
  //get_mobile()->load_positions(config.particles());
  DEBUG("selected " << mobile_.str());
  remove_unphysical_sites(config);
  ASSERT(mobile_.num_particles() > 0, "all sites shouldn't be unphysical");
  set_mobile_original(system);
  DEBUG("selected " << mobile_.str());
  return true;
}

std::shared_ptr<TrialSelect> TrialSelectAll::create(std::istream& istr) const {
  return std::make_shared<TrialSelectAll>(istr);
}

TrialSelectAll::TrialSelectAll(std::istream& istr)
  : TrialSelect(istr) {
  // ASSERT(class_name_ == "TrialSelectAll", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(1759 == version, "mismatch version: " << version);
}

void TrialSelectAll::serialize_trial_select_all_(std::ostream& ostr) const {
  serialize_trial_select_(ostr);
  feasst_serialize_version(1759, ostr);
}

void TrialSelectAll::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_select_all_(ostr);
}

}  // namespace feasst
