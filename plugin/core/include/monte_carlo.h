
#ifndef FEASST_CORE_MONTE_CARLO_H_
#define FEASST_CORE_MONTE_CARLO_H_

#include <memory>
#include "core/include/trial_factory.h"

namespace feasst {

class MonteCarlo {
 public:
  void set_criteria(std::shared_ptr<Criteria> criteria) { criteria_ = criteria; }
  void set_system(const System& system) { system_ = system; }
  void add_trial(std::shared_ptr<Trial> trial) {
    trial_factory_.add(trial);
  }
  const System& system() const { return system_; }
//  System * get_system() { return &system_; }
  std::shared_ptr<Criteria> get_criteria() { return criteria_; }
  const TrialFactory& trials() const { return trial_factory_; }

  void seek_num_particles(const int num) {
    TrialTransfer add;
    add.set_add_probability(1);
    while (system_.configuration().num_particles() < num) {
      attempt();
      add.attempt(criteria_.get(), &system_);
    }
  }

  void attempt(const int num_trials = 1) {
    for (int trial = 0; trial < num_trials; ++trial) {
      trial_factory_.attempt(criteria_.get(), &system_);
    }
  }

 private:
  std::shared_ptr<Criteria> criteria_;
  TrialFactory trial_factory_;
  System system_;
};

}  // namespace feasst

#endif  // FEASST_CORE_MONTE_CARLO_H_
