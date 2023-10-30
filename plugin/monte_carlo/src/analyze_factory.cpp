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

int AnalyzeFactory::min_block_(const Criteria& criteria) const {
  int min = 1e9;
  for (int state = criteria.soft_min(); state <= criteria.soft_max(); ++state) {
    const int num = static_cast<int>(
      analyzers_[state]->accumulator().blocks()[0].size());
    if (num < min) {
      min = num;
    }
  }
  return min;
}

std::string AnalyzeFactory::write_blocks_(const int min_block,
  const Accumulator& acc) const {
  std::stringstream ss;
  for (int block = 0; block < min_block; ++block) {
    const std::vector<std::vector<double>>& blocks = acc.blocks();
    if (static_cast<int>(blocks.size()) > 0) {
      if (static_cast<int>(blocks[0].size()) > block) {
        ss << blocks[0][block] << ",";
      } else {
        ss << ",";
      }
    } else {
      ss << ",";
    }
  }
  return ss.str();
}

void AnalyzeFactory::trial(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  DEBUG(" class? " << class_name());
  int stt = -1;
  if (is_multistate()) {
    stt = criteria.state();
  }
  DEBUG("criteria.num_iterations(" << stt << ") " << criteria.num_iterations(stt)  <<
       " start_after_iteration() " << start_after_iteration());
  if ( (stop_after_phase() != -1 && criteria.phase() > stop_after_phase()) ||
       (stop_after_iteration() != -1 && criteria.num_iterations(stt) > stop_after_iteration()) ||
       (criteria.phase() <= start_after_phase()) ||
       (criteria.num_iterations(stt) <= start_after_iteration()) ) {
    return;
  }
  if (is_multistate()) {
    DEBUG("multistate");
    DEBUG("state? " << criteria.state());
    if (is_multistate_aggregate()) {
      DEBUG("aggregating");
      DEBUG("sz " << analyzers_.size());
      DEBUG("state " << criteria.state());
      ASSERT(criteria.state() < static_cast<int>(analyzers_.size()),
        "state: " << criteria.state() << " >= multistate analyzers: " <<
        analyzers_.size() << ". Was a flat histogram simulation reinitialized"
        << " after a multistate Analyzer?");
      ASSERT(criteria.state() >= 0, "No state");
      analyzers_[criteria.state()]->check_update_(criteria, system, trial_factory);
      if (is_time(trials_per_write(), &trials_since_write_)) {
        const int min_block = min_block_(criteria);
        std::stringstream ss;
        for (int state = 0; state < num(); ++state) {
          if (state == 0) {
            ss << "state,";
            ss << analyzers_[state]->header(criteria, system, trial_factory);
            ss.seekp(-1, ss.cur); // remove endl
            for (int block = 0; block < min_block; ++block) {
              ss << "block" << block << ",";
            }
            ss << std::endl;
          }
          DEBUG("state " << state);
          DEBUG("crit " << criteria.state());
          DEBUG("crit " << criteria.num_states());
          ss << state << ",";
          ss << analyzers_[state]->write(criteria, system, trial_factory);
          ss.seekp(-1, ss.cur); // remove endl
          ss << write_blocks_(min_block, analyzers_[state]->accumulator());
          ss << std::endl;
        }
        printer(ss.str(), output_file(criteria));
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
      const int min_block = min_block_(criteria);
      std::stringstream ss;
      for (int state = 0; state < num(); ++state) {
        if (state == 0) {
          ss << "state,";
          ss << analyzers_[state]->header(criteria, system, trial_factory);
          ss.seekp(-1, ss.cur); // remove endl
          for (int block = 0; block < min_block; ++block) {
            ss << "block" << block << ",";
          }
          ss << std::endl;
        }
        ss << state << ",";
        ss << analyzers_[state]->write(criteria, system, trial_factory);
        ss.seekp(-1, ss.cur); // remove endl
        ss << write_blocks_(min_block, analyzers_[state]->accumulator());
        ss << std::endl;
      }
      printer(ss.str(), output_file(criteria));
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
