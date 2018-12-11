
#ifndef FEASST_CORE_TRIAL_FACTORY_H_
#define FEASST_CORE_TRIAL_FACTORY_H_

#include <memory>
#include "core/include/trial.h"

namespace feasst {

class TrialFactory : public Trial {
 public:
  virtual void attempt(Criteria* criteria, System * system) {
    increment_num_attempts();
    ASSERT(num_trials() > 0, "size error");
    const double random_uniform = random_.uniform();
    int attempt = 0;
    int index = 0;
    while (attempt == 0 && index < num_trials()) {
      if (random_uniform < cumulative_probability_[index]) {
        trials_[index]->attempt(criteria, system);
        trials_[index]->increment_num_attempts();
        attempt = 1;
        TRACE("attempt " << attempt);
      }
      ++index;
      TRACE("index " << index << " " << random_uniform);
    }
    ASSERT(attempt != 0, "trial should have been attempted");
  }

  void add(std::shared_ptr<Trial> trial) {
    trials_.push_back(trial);

    // update probability of selection
    double total_weight = 0.;
    for (std::shared_ptr<Trial> trial : trials_) {
      total_weight += trial->weight();
    }
    cumulative_probability_.clear();
    double probability = 0;
    for (std::shared_ptr<Trial> trial : trials_) {
      probability += trial->weight()/total_weight;
      cumulative_probability_.push_back(probability);
//      cout << cumulative_probability_.back() << endl;
    }
    ASSERT(std::abs(cumulative_probability_.back() - 1) < 1e-15, "weight error");
  }

  int num_trials() const { return trials_.size(); }

 private:
  std::vector<std::shared_ptr<Trial> > trials_;
  std::vector<double> cumulative_probability_;
  Random random_;
};

}  // namespace feasst

#endif  // FEASST_CORE_TRIAL_FACTORY_H_
