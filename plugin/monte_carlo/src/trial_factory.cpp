#include "monte_carlo/include/trial_factory.h"

namespace feasst {

TrialFactory::TrialFactory() { class_name_ = "TrialFactory"; }

class MapTrialFactory {
 public:
  MapTrialFactory() {
    auto obj = MakeTrialFactory();
    obj->deserialize_map()["TrialFactory"] = obj;
  }
};

static MapTrialFactory mapper_ = MapTrialFactory();

std::shared_ptr<Trial> TrialFactory::create(std::istream& istr) const {
  return std::make_shared<TrialFactory>(istr);
}

TrialFactory::TrialFactory(std::istream& istr) : Trial(istr) {
  // ASSERT(class_name_ == "TrialFactory", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(189 == version, "mismatch version: " << version);
  // HWH for unknown reasons, this function template does not work.
  // feasst_deserialize(&trials_, istr);
  { int dim1;
    istr >> dim1;
    trials_.resize(dim1);
    for (int index = 0; index < dim1; ++index) {
      // feasst_deserialize_fstobj((trials_)[index], istr);
      int existing;
      istr >> existing;
      if (existing != 0) {
        trials_[index] = trials_[index]->deserialize(istr);
      }
    }
  }
  feasst_deserialize(&cumulative_probability_, istr);
}

void TrialFactory::serialize_trial_factory_(std::ostream& ostr) const {
  serialize_trial_(ostr);
  feasst_serialize_version(189, ostr);
  feasst_serialize(trials_, ostr);
  feasst_serialize(cumulative_probability_, ostr);
}

void TrialFactory::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_factory_(ostr);
}

}  // namespace feasst
