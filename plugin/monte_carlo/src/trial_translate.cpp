#include "monte_carlo/include/trial_translate.h"

namespace feasst {

class MapTrialTranslate {
 public:
  MapTrialTranslate() {
    auto obj = MakeTrialTranslate();
    obj->deserialize_map()["TrialTranslate"] = obj;
  }
};

static MapTrialTranslate mapper_ = MapTrialTranslate();

std::shared_ptr<Trial> TrialTranslate::create(std::istream& istr) const {
  return std::make_shared<TrialTranslate>(istr);
}

TrialTranslate::TrialTranslate(std::istream& istr) : TrialMove(istr) {
  // ASSERT(class_name_ == "TrialTranslate", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(266 == version, "mismatch version: " << version);
}

void TrialTranslate::serialize_trial_translate_(std::ostream& ostr) const {
  serialize_trial_move_(ostr);
  feasst_serialize_version(266, ostr);
}

void TrialTranslate::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_translate_(ostr);
}

}  // namespace feasst
