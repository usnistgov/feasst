
#ifndef FEASST_MONTE_CARLO_TRIAL_FACTORY_H_
#define FEASST_MONTE_CARLO_TRIAL_FACTORY_H_

#include <sstream>
#include <memory>
#include "monte_carlo/include/trial.h"
// #include "utils/include/timer.h"

namespace feasst {

/**
  Contains multiple Trials.
 */
class TrialFactory : public Trial {
 public:
  TrialFactory();

  /// Add a trial.
  void add(std::shared_ptr<Trial> trial);

  /// Return the number of trials.
  int num() const { return static_cast<int>(trials_.size()); }

  /// Return a trial by index of the order trials were added.
  const Trial& trial(const int index) const override {
    return const_cast<Trial&>(*trials_[index]); }

  const std::vector<std::shared_ptr<Trial> >& trials() const override {
    return trials_; }

  /// Return the index of a trial selected with probability proportional to its
  /// weight.
  int random_index(Random * random);

  /// Attempt one of the trials. Return true if accepted.
  bool attempt(
      Criteria* criteria,
      System * system,
      /// attempt trial_index. If -1, choose randomly with probabilty
      /// determined from the weight.
      int trial_index,
      Random * random);

  /// Attempt one of the trials with selection probability proportional to
  /// the weight.
  bool attempt(Criteria* criteria, System * system, Random * random) override {
    return attempt(criteria, system, -1, random); }

  /// Revert changes to system by trial index.
  void revert(const int index, const bool accepted, System * system);

  void finalize(const int index, System * system) {
    trials_[index]->finalize(system); }

  /// Return the header description for the statuses of the trials (e.g.,
  /// acceptance, etc).
  std::string status_header() const override;

  // HWH hackish interface for prefetch
  // Require manual finalization of trials (e.g., Prefetch).
  void delay_finalize();
  void imitate_trial_rejection_(const int index);

  /// Return the statuses of the trials (e.g., acceptance, etc).
  std::string status() const override;

  void reset_stats() override;
  void tune() override;
  void precompute(Criteria * criteria, System * system) override;

//  const Timer& timer() const { return timer_; }

  bool is_equal(const TrialFactory& factory) const;

  // serialize
  std::shared_ptr<Trial> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialFactory(std::istream& istr);
  virtual ~TrialFactory() {}

 protected:
  void serialize_trial_factory_(std::ostream& ostr) const;

 private:
  std::vector<std::shared_ptr<Trial> > trials_;
  std::vector<double> cumulative_probability_;
//  Timer timer_;
};

inline std::shared_ptr<TrialFactory> MakeTrialFactory() {
  return std::make_shared<TrialFactory>();
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_FACTORY_H_
