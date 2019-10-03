#include "monte_carlo/include/trial_add.h"

namespace feasst {

TrialAdd::TrialAdd(const argtype& args) : Trial(args) {
  add_stage(
    std::make_shared<TrialSelectParticle>(args),
    std::make_shared<PerturbAdd>(args),
    args
  );
  set(std::make_shared<TrialComputeAdd>());
  class_name_ = "TrialAdd";
}

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
