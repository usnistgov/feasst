#include "utils/include/serialize.h"
#include "utils/include/timer_rdtsc.h"
#include "math/include/accumulator.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_factory.h"
#include "monte_carlo/include/monte_carlo.h"
#include "monte_carlo/include/analyze_factory.h"

namespace feasst {

FEASST_MAPPER(AnalyzeFactory,);

AnalyzeFactory::AnalyzeFactory(argtype args) : Analyze(&args) {}
AnalyzeFactory::~AnalyzeFactory() {}

void AnalyzeFactory::initialize(MonteCarlo * mc) {
  for (std::shared_ptr<Analyze> analyze : analyzers_) {
    analyze->initialize(mc);
  }
}

void AnalyzeFactory::add(std::shared_ptr<Analyze> analyze) {
  analyzers_.push_back(analyze);
  if (timer_) timer_->add();
}

void AnalyzeFactory::remove(const int index) {
  analyzers_.erase(analyzers_.begin() + index);
  if (timer_) timer_->erase(index);
}

void AnalyzeFactory::trial_(const MonteCarlo& mc, const int index) {
  DEBUG("index " << index << " sz " << analyzers_.size());
  ASSERT(index < static_cast<int>(analyzers_.size()),
    "index: " << index << " too large when there are " << analyzers_.size());
  DEBUG(analyzers_[index]->class_name());
  if (timer_) timer_->start(index);
  analyzers_[index]->trial(mc);
  if (timer_) timer_->start(-1);
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

void AnalyzeFactory::trial(const MonteCarlo& mc) {
  DEBUG(" class? " << class_name());
  const Criteria& criteria = mc.criteria();
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
      analyzers_[criteria.state()]->check_update_(mc);
      if (is_time(trials_per_write(), &trials_since_write_)) {
        const int min_block = min_block_(criteria);
        std::stringstream ss;
        for (int state = 0; state < num(); ++state) {
          const Accumulator acc = analyzers_[state]->accumulator();
          const bool acc_used = acc.num_values() != 0;
          if (state == 0) {
            ss << "state,";
            ss << analyzers_[state]->header(mc);
            if (acc_used) {
              ss.seekp(-1, ss.cur); // remove endl
              for (int block = 0; block < min_block; ++block) {
                ss << "block" << block << ",";
              }
              ss << std::endl;
            }
          }
          DEBUG("state " << state);
          DEBUG("crit " << criteria.state());
          DEBUG("crit " << criteria.num_states());
          ss << state << ",";
          ss << analyzers_[state]->write(mc);
          if (acc_used) {
            ss.seekp(-1, ss.cur); // remove endl
            ss << write_blocks_(min_block, analyzers_[state]->accumulator());
            ss << std::endl;
          }
        }
        printer(ss.str(), output_file(criteria));
      }
    } else {
      trial_(mc, criteria.state());
    }
  } else {
    DEBUG("not multistate");
    for (int index = 0; index < num(); ++index) {
//    for (const std::shared_ptr<Analyze> analyze : analyzers_) {
      DEBUG("index " << index);
      trial_(mc, index);
    }
  }
}

void AnalyzeFactory::write_to_file(const MonteCarlo& mc) {
  const Criteria& criteria = mc.criteria();
  if (is_multistate()) {
    if (is_multistate_aggregate()) {
      const int min_block = min_block_(criteria);
      std::stringstream ss;
      for (int state = 0; state < num(); ++state) {
        if (state == 0) {
          ss << "state,";
          ss << analyzers_[state]->header(mc);
          ss.seekp(-1, ss.cur); // remove endl
          for (int block = 0; block < min_block; ++block) {
            ss << "block" << block << ",";
          }
          ss << std::endl;
        }
        ss << state << ",";
        ss << analyzers_[state]->write(mc);
        ss.seekp(-1, ss.cur); // remove endl
        ss << write_blocks_(min_block, analyzers_[state]->accumulator());
        ss << std::endl;
      }
      printer(ss.str(), output_file(criteria));
    } else {
      analyzers_[criteria.state()]->write_to_file(mc);
    }
  } else {
    for (int index = 0; index < num(); ++index) {
      analyzers_[index]->write_to_file(mc);
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

void AnalyzeFactory::set_timer() {
  timer_ = std::make_unique<TimerRDTSC>(num());
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
  feasst_deserialize(timer_, istr);
}

void AnalyzeFactory::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(1640, ostr);
  feasst_serialize_fstdr(analyzers_, ostr);
  feasst_serialize(timer_, ostr);
}

}  // namespace feasst
