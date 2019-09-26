#include "chain/include/trial_select_end_segment.h"

namespace feasst {

class MapTrialSelectEndSegment {
 public:
  MapTrialSelectEndSegment() {
    auto obj = MakeTrialSelectEndSegment();
    obj->deserialize_map()["TrialSelectEndSegment"] = obj;
  }
};

static MapTrialSelectEndSegment mapper_ = MapTrialSelectEndSegment();

std::shared_ptr<TrialSelect> TrialSelectEndSegment::create(std::istream& istr) const {
  return std::make_shared<TrialSelectEndSegment>(istr);
}

TrialSelectEndSegment::TrialSelectEndSegment(std::istream& istr)
  : TrialSelectSegment(istr) {
  // ASSERT(class_name_ == "TrialSelectEndSegment", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(456 == version, "mismatch version: " << version);
}

void TrialSelectEndSegment::serialize_trial_select_end_segment_(std::ostream& ostr) const {
  serialize_trial_select_segment_(ostr);
  feasst_serialize_version(456, ostr);
}

void TrialSelectEndSegment::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_select_end_segment_(ostr);
}

}  // namespace feasst
