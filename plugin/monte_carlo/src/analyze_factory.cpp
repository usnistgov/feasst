#include <iostream>
#include "monte_carlo/include/analyze_factory.h"
#include "utils/include/serialize.h"

namespace feasst {

class MapAnalyzeFactory {
 public:
  MapAnalyzeFactory() {
    AnalyzeFactory().deserialize_map()["AnalyzeFactory"] =
      std::make_shared<AnalyzeFactory>();
  }
};

static MapAnalyzeFactory mapper_ = MapAnalyzeFactory();

void AnalyzeFactory::initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  for (std::shared_ptr<Analyze> analyze : analyzers_) {
    analyze->initialize(criteria, system, trial_factory);
  }
}

void AnalyzeFactory::trial(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  DEBUG(" class? " << class_name());
  if (is_multistate()) {
    DEBUG("multistate");
    DEBUG("state? " << criteria.state());
    if (is_multistate_aggregate()) {
      DEBUG("aggregating");
      analyzers_[criteria.state()]->check_update_(criteria, system, trial_factory);
      if (is_time(steps_per_write(), &steps_since_write_)) {
        std::stringstream ss;
        for (int state = 0; state < num(); ++state) {
          if (state == 0) {
            ss << "state,";
            ss << analyzers_[state]->header(criteria, system, trial_factory);
          }
          DEBUG("state " << state);
          DEBUG("crit " << criteria.state());
          DEBUG("crit " << criteria.num_states());
          ss << state << ",";
          ss << analyzers_[state]->write(criteria, system, trial_factory);
//          ss << std::endl;
        }
        printer(ss.str());
      }
    } else {
      trial_(criteria, system, trial_factory, criteria.state());
    }
  } else {
    DEBUG("not multistate");
    for (int index = 0; index < num(); ++index) {
//    for (const std::shared_ptr<Analyze> analyze : analyzers_) {
      DEBUG("index " << index);
      trial_(criteria, system, trial_factory, index);
    }
  }
}

void AnalyzeFactory::trial_(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory,
    const int index) {
  // timer_.start(index + 1);
  DEBUG("index " << index << " sz " << analyzers_.size());
  ASSERT(index < static_cast<int>(analyzers_.size()),
    "index: " << index << " too large when there are " << analyzers_.size());
  DEBUG(analyzers_[index]->class_name());
  analyzers_[index]->trial(criteria, system, trial_factory);
  // timer_.end();
}

AnalyzeFactory::AnalyzeFactory(std::istream& istr) : Analyze(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 1640, "unrecognized verison: " << version);
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
  feasst_serialize_version(1640, ostr);
  feasst_serialize_fstdr(analyzers_, ostr);
}

}  // namespace feasst
