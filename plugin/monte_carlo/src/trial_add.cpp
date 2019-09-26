#include "monte_carlo/include/trial_add.h"

namespace feasst {

class MapTrialAdd {
 public:
  MapTrialAdd() {
    auto obj = MakeTrialAdd();
    obj->deserialize_map()["TrialAdd"] = obj;
  }
};

static MapTrialAdd mapper_ = MapTrialAdd();

std::shared_ptr<Trial> TrialAdd::create(std::istream& istr) const {
  return std::make_shared<TrialAdd>(istr);
}

TrialAdd::TrialAdd(std::istream& istr) : Trial(istr) {
  // ASSERT(class_name_ == "TrialAdd", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(736 == version, "mismatch version: " << version);
}

void TrialAdd::serialize_trial_add_(std::ostream& ostr) const {
  serialize_trial_(ostr);
  feasst_serialize_version(736, ostr);
}

void TrialAdd::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_add_(ostr);
}

}  // namespace feasst
