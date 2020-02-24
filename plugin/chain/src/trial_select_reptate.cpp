#include "chain/include/trial_select_reptate.h"

namespace feasst {

class MapTrialSelectReptate {
 public:
  MapTrialSelectReptate() {
    auto obj = MakeTrialSelectReptate({{"max_length", "1"}});
    obj->deserialize_map()["TrialSelectReptate"] = obj;
  }
};

static MapTrialSelectReptate mapper_ = MapTrialSelectReptate();

std::shared_ptr<TrialSelect> TrialSelectReptate::create(std::istream& istr) const {
  return std::make_shared<TrialSelectReptate>(istr);
}

TrialSelectReptate::TrialSelectReptate(std::istream& istr)
  : TrialSelectEndSegment(istr) {
  // ASSERT(class_name_ == "TrialSelectReptate", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(812 == version, "mismatch version: " << version);
}

void TrialSelectReptate::serialize_trial_select_reptate_(std::ostream& ostr) const {
  serialize_trial_select_end_segment_(ostr);
  feasst_serialize_version(812, ostr);
}

void TrialSelectReptate::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_select_reptate_(ostr);
}

}  // namespace feasst
