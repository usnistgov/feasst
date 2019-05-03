
#include "core/include/analyze_factory.h"

namespace feasst {

class MapAnalyzeFactory {
 public:
  MapAnalyzeFactory() {
    AnalyzeFactory().deserialize_map()["AnalyzeFactory"] =
      std::make_shared<AnalyzeFactory>();
  }
};

static MapAnalyzeFactory mapper_ = MapAnalyzeFactory();

void AnalyzeFactory::initialize(const std::shared_ptr<Criteria> criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  for (std::shared_ptr<Analyze> analyze : analyzers_) {
    analyze->initialize(criteria, system, trial_factory);
  }
}

void AnalyzeFactory::trial(const std::shared_ptr<Criteria> criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  DEBUG("multistate? " << is_multistate() << " class? " << class_name());
  if (is_multistate()) {
    DEBUG("state? " << criteria->state() << " class " << analyzers_[criteria->state()]->class_name());
    analyzers_[criteria->state()]->trial(criteria, system, trial_factory);
  } else {
    for (const std::shared_ptr<Analyze> analyze : analyzers_) {
      DEBUG(analyze->class_name());
      analyze->trial(criteria, system, trial_factory);
    }
  }
}

AnalyzeFactory::AnalyzeFactory(std::istream& istr) : Analyze(istr) {
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
  Stepper::serialize(ostr);
  feasst_serialize_version(1, ostr);
  feasst_serialize_fstdr(analyzers_, ostr);
}

}  // namespace feasst
