
#ifndef FEASST_CORE_TRIAL_FACTORY_H_
#define FEASST_CORE_TRIAL_FACTORY_H_

#include <memory>
#include "core/include/trial.h"

namespace feasst {

class TrialFactory : public Trial {
 public:
  void attempt(
      Criteria* criteria,
      System * system) override {
    attempt(criteria, system, -1);
  }
  void attempt(
      Criteria* criteria,
      System * system,
      /// attempt trial_index. If -1, choose randomly with probabilty
      /// determined from the weight.
      const int trial_index) {
    increment_num_attempts();
    ASSERT(num_trials() > 0, "size error");
    if (trial_index != -1) {
      attempt_(criteria, system, trial_index);
      return;
    }
    const double random_uniform = random_.uniform();
    int attempt = 0;
    int index = 0;
    while (attempt == 0 && index < num_trials()) {
      if (random_uniform < cumulative_probability_[index]) {
        attempt_(criteria, system, index);
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

  std::vector<std::shared_ptr<Trial> > trials() { return trials_; }

  /// Return the header description for the statuses of the trials (e.g., acceptance, etc).
  std::string status_header() const override {
    std::stringstream ss;
    ss << "attempt ";
    for (const std::shared_ptr<Trial> trial : trials_) {
      ss << trial->status_header() << " ";
    }
    return ss.str();
  }

  /// Return the statuses of the trials (e.g., acceptance, etc).
  std::string status() const override {
    std::stringstream ss;
    ss << num_attempts() << " ";
    for (const std::shared_ptr<Trial> trial : trials_) {
      ss << trial->status() << " ";
    }
    return ss.str();
  }

  void reset_stats() override {
    Trial::reset_stats();
    for (std::shared_ptr<Trial> trial : trials_) {
      trial->reset_stats();
    }
  }

  void tune() override {
    for (std::shared_ptr<Trial> trial : trials_) {
      trial->tune();
    }
  }

 private:
  std::vector<std::shared_ptr<Trial> > trials_;
  std::vector<double> cumulative_probability_;
  Random random_;

  void attempt_(
      Criteria* criteria,
      System * system,
      const int index) {
    trials_[index]->attempt(criteria, system);
  }
};

}  // namespace feasst

#endif  // FEASST_CORE_TRIAL_FACTORY_H_
