#include "chain/include/trial_select_segment.h"

namespace feasst {

class MapTrialSelectSegment {
 public:
  MapTrialSelectSegment() {
    auto obj = MakeTrialSelectSegment();
    obj->deserialize_map()["TrialSelectSegment"] = obj;
  }
};

static MapTrialSelectSegment mapper_ = MapTrialSelectSegment();

std::shared_ptr<TrialSelect> TrialSelectSegment::create(std::istream& istr) const {
  return std::make_shared<TrialSelectSegment>(istr);
}

TrialSelectSegment::TrialSelectSegment(std::istream& istr)
  : TrialSelectParticle(istr) {
  // ASSERT(class_name_ == "TrialSelectSegment", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(658 == version, "mismatch version: " << version);
  feasst_deserialize(&max_length_, istr);
}

void TrialSelectSegment::serialize_trial_select_segment_(std::ostream& ostr) const {
  serialize_trial_select_particle_(ostr);
  feasst_serialize_version(658, ostr);
  feasst_serialize(max_length_, ostr);
}

void TrialSelectSegment::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_select_segment_(ostr);
}

}  // namespace feasst
