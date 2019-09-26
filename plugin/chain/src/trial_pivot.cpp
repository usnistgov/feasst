#include "chain/include/trial_pivot.h"

namespace feasst {

class MapTrialPivot {
 public:
  MapTrialPivot() {
    auto obj = MakeTrialPivot();
    obj->deserialize_map()["TrialPivot"] = obj;
  }
};

static MapTrialPivot mapper_ = MapTrialPivot();

std::shared_ptr<Trial> TrialPivot::create(std::istream& istr) const {
  return std::make_shared<TrialPivot>(istr);
}

TrialPivot::TrialPivot(std::istream& istr) : TrialMove(istr) {
  // ASSERT(class_name_ == "TrialPivot", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(718 == version, "mismatch version: " << version);
}

void TrialPivot::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_move_(ostr);
  feasst_serialize_version(718, ostr);
}

}  // namespace feasst
