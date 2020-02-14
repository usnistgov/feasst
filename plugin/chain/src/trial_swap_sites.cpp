#include "chain/include/trial_swap_sites.h"

namespace feasst {

class MapTrialSwapSites {
 public:
  MapTrialSwapSites() {
    auto obj = MakeTrialSwapSites({{"site_type1", "0"}, {"site_type2", "1"}});
    obj->deserialize_map()["TrialSwapSites"] = obj;
  }
};

static MapTrialSwapSites mapper_ = MapTrialSwapSites();

std::shared_ptr<Trial> TrialSwapSites::create(std::istream& istr) const {
  return std::make_shared<TrialSwapSites>(istr);
}

TrialSwapSites::TrialSwapSites(std::istream& istr) : Trial(istr) {
  // ASSERT(class_name_ == "TrialSwapSites", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(346 == version, "mismatch version: " << version);
}

void TrialSwapSites::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_(ostr);
  feasst_serialize_version(346, ostr);
}

}  // namespace feasst
