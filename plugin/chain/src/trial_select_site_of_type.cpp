#include "chain/include/trial_select_site_of_type.h"

namespace feasst {

class MapTrialSelectSiteOfType {
 public:
  MapTrialSelectSiteOfType() {
    auto obj = MakeTrialSelectSiteOfType({{"site_type", "0"}});
    obj->deserialize_map()["TrialSelectSiteOfType"] = obj;
  }
};

static MapTrialSelectSiteOfType mapper_ = MapTrialSelectSiteOfType();

std::shared_ptr<TrialSelect> TrialSelectSiteOfType::create(std::istream& istr) const {
  return std::make_shared<TrialSelectSiteOfType>(istr);
}

TrialSelectSiteOfType::TrialSelectSiteOfType(std::istream& istr)
  : TrialSelect(istr) {
  // ASSERT(class_name_ == "TrialSelectSiteOfType", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(272 == version, "mismatch version: " << version);
  feasst_deserialize(&site_type_, istr);
}

void TrialSelectSiteOfType::serialize_trial_select_segment_(std::ostream& ostr) const {
  serialize_trial_select_(ostr);
  feasst_serialize_version(272, ostr);
  feasst_serialize(site_type_, ostr);
}

void TrialSelectSiteOfType::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_select_segment_(ostr);
}

}  // namespace feasst
