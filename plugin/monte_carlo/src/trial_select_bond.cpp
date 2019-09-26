#include "monte_carlo/include/trial_select_bond.h"

namespace feasst {

class MapTrialSelectBond {
 public:
  MapTrialSelectBond() {
    auto obj = MakeTrialSelectBond({{"mobile_site", "1"}, {"anchor_site", "0"}});
    obj->deserialize_map()["TrialSelectBond"] = obj;
  }
};

static MapTrialSelectBond mapper_ = MapTrialSelectBond();

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
