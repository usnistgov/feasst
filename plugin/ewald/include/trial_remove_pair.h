
#ifndef FEASST_MONTE_CARLO_TRIAL_REMOVE_PAIR_H_
#define FEASST_MONTE_CARLO_TRIAL_REMOVE_PAIR_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"

namespace feasst {

/// Attempt to remove a particle.
class TrialRemovePair : public Trial {
 public:
  explicit TrialRemovePair(const argtype& args = argtype());
  std::shared_ptr<Trial> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialRemovePair(std::istream& istr);
  virtual ~TrialRemovePair() {}

 protected:
  void serialize_trial_remove_pair_(std::ostream& ostr) const;
};

inline std::shared_ptr<TrialRemovePair> MakeTrialRemovePair(
    const argtype &args = argtype()) {
  return std::make_shared<TrialRemovePair>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_REMOVE_PAIR_H_
