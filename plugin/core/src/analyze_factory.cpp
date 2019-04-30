
#include "core/include/analyze_factory.h"

namespace feasst {

void AnalyzeFactory::initialize(const std::shared_ptr<Criteria> criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  for (const std::shared_ptr<Analyze> analyze : analyzers_) {
    analyze->initialize(criteria, system, trial_factory);
  }
}

void AnalyzeFactory::trial(const std::shared_ptr<Criteria> criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  for (const std::shared_ptr<Analyze> analyze : analyzers_) {
    analyze->trial(criteria, system, trial_factory);
  }
}

AnalyzeFactory::AnalyzeFactory(std::istream& istr) {
  std::string class_name;
  istr >> class_name;
  ASSERT(class_name == class_name_, "class_name(" << class_name << ") does "
    << "not match class_name_(" << class_name_ << ")");
  feasst_deserialize_version(istr);
  // feasst_deserialize_fstdr(&analyzers_, istr);
  // HWH for unknown reasons, function template doesn't work
  int dim1;
  istr >> dim1;
  analyzers_.resize(dim1);
  for (int index = 0; index < dim1; ++index) {
    int existing;
    istr >> existing;
    if (existing != 0) {
      analyzers_[index] = analyzers_[index]->deserialize(istr);
    }
  }
}

void AnalyzeFactory::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  feasst_serialize_version(1, ostr);
  feasst_serialize_fstdr(analyzers_, ostr);
}

}  // namespace feasst
