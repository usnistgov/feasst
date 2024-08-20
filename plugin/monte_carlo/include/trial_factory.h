
#ifndef FEASST_MONTE_CARLO_TRIAL_FACTORY_H_
#define FEASST_MONTE_CARLO_TRIAL_FACTORY_H_

#include <map>
#include <memory>
#include <string>
#include <vector>
#include "monte_carlo/include/trial.h"
// #include "utils/include/timer.h"

namespace feasst {

class Random;

typedef std::map<std::string, std::string> argtype;

/**
  Contains multiple Trials.
 */
class TrialFactory : public Trial {
 public:
  explicit TrialFactory(argtype args = argtype());
  explicit TrialFactory(argtype * args);

  /// Add a trial.
  void add(std::shared_ptr<Trial> trial);

  /// Remove a trial by index.
  void remove(const int index);

  /// Return the number of trials.
  int num() const { return static_cast<int>(trials_.size()); }

  /// Return a trial by index of the order trials were added.
  const Trial& trial(const int index) const override;

  const std::vector<std::shared_ptr<Trial> >& trials() const override {
    return trials_; }

  /// Return the index of a trial selected with probability proportional to its
  /// weight.
  int random_index(Random * random);

  /// Attempt one of the trials. Return true if accepted.
  bool attempt(
    Criteria * criteria,
    System * system,
    /// attempt trial_index. If -1, choose randomly with probabilty
    /// determined from the weight.
    int trial_index,
    Random * random);

  /// Attempt one of the trials with selection probability proportional to
  /// the weight.
  bool attempt(Criteria * criteria, System * system, Random * random) override {
    return attempt(criteria, system, -1, random); }

  /// Return the index of the last trial attempted.
  int last_index() const { return last_index_; }
  void set_last_index(const int index) { last_index_ = index; }

  /// Return the cumulative probability of each trial.
  const std::vector<double>& cumulative_probability() const {
    return data_.dble_2D()[0]; }

  /// Revert changes to system by trial index.
  void revert(const int index, const bool accepted, const bool auto_rejected,
    System * system, Criteria * criteria);

  void finalize(const int index, System * system, Criteria * criteria) {
    trials_[index]->finalize(system, criteria); }

  /// Return the header description for the statuses of the trials (e.g.,
  /// acceptance, etc).
  std::string status_header() const override;

  // HWH hackish interface for prefetch
  // Require manual finalization of trials (e.g., Prefetch).
  void delay_finalize();
  void imitate_trial_rejection_(const int index,
    const bool auto_reject);

  /// Return the statuses of the trials (e.g., acceptance, etc).
  std::string status() const override;

  void reset_stats() override;
  void tune() override;
  void precompute(Criteria * criteria, System * system) override;

//  const Timer& timer() const { return timer_; }

  bool is_equal(const TrialFactory& factory) const;
  void synchronize_(const Trial& trial) override;
  Trial * get_trial(const int index) { return trials_[index].get(); }

  void set_tunable(const int trial_index, const double tunable);

  // serialize
  std::shared_ptr<Trial> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialFactory(std::istream& istr);
  virtual ~TrialFactory() {}

 protected:
  void serialize_trial_factory_(std::ostream& ostr) const;

 private:
  std::vector<std::shared_ptr<Trial> > trials_;
  std::vector<double> * get_cumulative_probability_() {
    return &((*data_.get_dble_2D())[0]); }
  // std::vector<double> cumulative_probability_;
  bool adjustable_weights_ = false;

  // not to be serialized
  int last_index_ = -1;
//  Timer timer_;

  void update_cumul_prob_();
};

inline std::shared_ptr<TrialFactory> MakeTrialFactory() {
  return std::make_shared<TrialFactory>();
}

/**
  Contains multiple Trials for use as input into MonteCarlo.
  As opposed to the above, these classes are named to enable factories.
  Serialization is not required.
 */
class TrialFactoryNamed {
 public:
  TrialFactoryNamed() {}
  const std::vector<std::shared_ptr<Trial> >& trials() const { return trials_; }
  // serialize
  std::string class_name() const { return class_name_; }
  void add(std::shared_ptr<Trial> trial) { trials_.push_back(trial); }
  void precompute(Criteria * criteria, System * system);
  Trial * trial_(int index) { return trials_[index].get(); }
//  virtual void serialize(std::ostream& ostr) const;
//  virtual std::shared_ptr<TrialFactoryNamed> create(std::istream& istr) const;
  virtual std::shared_ptr<TrialFactoryNamed> create(argtype * args) const;
  std::map<std::string, std::shared_ptr<TrialFactoryNamed> >& deserialize_map();
//  std::shared_ptr<TrialFactoryNamed> deserialize(std::istream& istr);
  std::shared_ptr<TrialFactoryNamed> factory(
    const std::string name, argtype * args);
//  explicit TrialFactoryNamed(std::istream& istr);
  virtual ~TrialFactoryNamed() {}

 protected:
  std::string class_name_;

 private:
  std::vector<std::shared_ptr<Trial> > trials_;
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_FACTORY_H_
