
#include "monte_carlo/include/modify_factory.h"

namespace feasst {

void ModifyFactory::initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  for (const std::shared_ptr<Modify> modify : modifiers_) {
    modify->initialize(criteria, system, trial_factory);
  }
}

void ModifyFactory::trial(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  //for (const std::shared_ptr<Modify> modify : modifiers_) {
  for (int index = 0; index < static_cast<int>(modifiers_.size()); ++index) {
    modifiers_[index]->trial(criteria, system, trial_factory);
  }
}

ModifyFactory::ModifyFactory(std::istream& istr) : Modify(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 7177, "unrecognized verison: " << version);
  // feasst_deserialize_fstdr(&modifiers_, istr);
  // HWH for unknown reasons, function template doesn't work
  int dim1;
  istr >> dim1;
  modifiers_.resize(dim1);
  for (int index = 0; index < dim1; ++index) {
    int existing;
    istr >> existing;
    if (existing != 0) {
      modifiers_[index] = modifiers_[index]->deserialize(istr);
    }
  }
}

void ModifyFactory::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(7177, ostr);
  feasst_serialize_fstdr(modifiers_, ostr);
}

}  // namespace feasst
