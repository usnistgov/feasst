
#ifndef FEASST_CORE_MONTE_CARLO_H_
#define FEASST_CORE_MONTE_CARLO_H_

#include <vector>
#include <memory>
#include "core/include/trial_factory.h"
#include "core/include/trial_transfer.h"
#include "core/include/analyze.h"
#include "core/include/modify.h"

namespace feasst {

class MonteCarlo {
 public:
  void set_criteria(std::shared_ptr<Criteria> criteria) { criteria_ = criteria; }
  void set_system(const System& system) { system_ = system; }
  void add_trial(std::shared_ptr<Trial> trial) {
    trial_factory_.add(trial);
  }
  const System& system() const { return system_; }
  System * get_system() { return &system_; }  // HWH depreciate: testing only
  std::shared_ptr<Criteria> get_criteria() { return criteria_; }
  const TrialFactory& trials() const { return trial_factory_; }

  void seek_num_particles(const int num) {
    TrialTransfer add;
    add.set_add_probability(1);
    while (system_.configuration().num_particles() < num) {
      attempt();
      add.attempt(criteria_.get(), &system_);
    }
    trial_factory_.reset_stats();
  }

  void attempt(const int num_trials = 1) {
    for (int trial = 0; trial < num_trials; ++trial) {
      trial_factory_.attempt(criteria_.get(), &system_);
      analyze_factory_.trial(criteria_, system_, trial_factory_);
      modify_factory_.trial(criteria_, &system_, &trial_factory_);
    }
  }

  void add_analyze(const std::shared_ptr<Analyze> analyze) {
    analyze->initialize(criteria_, system_, trial_factory_);
    analyze_factory_.add(analyze);
  }

  void add_modify(const std::shared_ptr<Modify> modify) {
    modify->initialize(criteria_, &system_, &trial_factory_);
    modify_factory_.add(modify);
  }

 private:
  std::shared_ptr<Criteria> criteria_;
  TrialFactory trial_factory_;
  System system_;
  AnalyzeFactory analyze_factory_;
  ModifyFactory modify_factory_;
};

}  // namespace feasst

#endif  // FEASST_CORE_MONTE_CARLO_H_
