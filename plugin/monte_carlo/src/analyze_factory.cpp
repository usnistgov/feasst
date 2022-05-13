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

void AnalyzeFactory::trial(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  DEBUG(" class? " << class_name());
  if ( (stop_after_phase() != -1 &&
        criteria.phase() > stop_after_phase()) ||
       (criteria.phase() <= start_after_phase()) ) {
    return;
  }
  if (is_multistate()) {
    DEBUG("multistate");
    DEBUG("state? " << criteria.state());
    if (is_multistate_aggregate()) {
      DEBUG("aggregating");
      DEBUG("sz " << analyzers_.size());
      ASSERT(criteria.state() < static_cast<int>(analyzers_.size()),
        "state: " << criteria.state() << " >= multistate analyzers: " <<
        analyzers_.size() << ". Was a flat histogram simulation reinitialized"
        << " after a multistate Analyzer?");
      analyzers_[criteria.state()]->check_update_(criteria, system, trial_factory);
      if (is_time(trials_per_write(), &trials_since_write_)) {
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
        printer(ss.str(), file_name(criteria));
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

void AnalyzeFactory::write_to_file(const Criteria& criteria,
  const System& system,
  const TrialFactory& trial_factory) {
  if (is_multistate()) {
    if (is_multistate_aggregate()) {
      std::stringstream ss;
      for (int state = 0; state < num(); ++state) {
        if (state == 0) {
          ss << "state,";
          ss << analyzers_[state]->header(criteria, system, trial_factory);
        }
        ss << state << ",";
        ss << analyzers_[state]->write(criteria, system, trial_factory);
      }
      printer(ss.str(), file_name(criteria));
    } else {
      analyzers_[criteria.state()]->write_to_file(criteria, system, trial_factory);
    }
  } else {
    for (int index = 0; index < num(); ++index) {
      analyzers_[index]->write_to_file(criteria, system, trial_factory);
    }
  }
}

void AnalyzeFactory::adjust_bounds(const bool adjusted_up,
    const std::vector<int>& states,
    AnalyzeFactory * analyze_factory) {
  for (int ai = 0; ai < num(); ++ai) {
    if (analyze(ai).is_multistate()) {
      DEBUG("ai " << ai);
      for (const int state : states) {
        if (adjusted_up) {
          *analyze_factory->get_analyze(ai)->get_analyze(state) = analyze(ai).analyze(state);
        } else {
          *get_analyze(ai)->get_analyze(state) = analyze_factory->analyze(ai).analyze(state);
        }
      }
    }
  }
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
