
#ifndef FEASST_MONTE_CARLO_TRIAL_REMOVE_MULTIPLE_H_
#define FEASST_MONTE_CARLO_TRIAL_REMOVE_MULTIPLE_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"

namespace feasst {

/// Attempt to remove a pair of particles.
class TrialRemoveMultiple : public Trial {
 public:
  explicit TrialRemoveMultiple(const argtype& args = argtype());
  std::shared_ptr<Trial> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialRemoveMultiple(std::istream& istr);
  virtual ~TrialRemoveMultiple() {}

 protected:
  void serialize_trial_remove_multiple_(std::ostream& ostr) const;
};

inline std::shared_ptr<TrialRemoveMultiple> MakeTrialRemoveMultiple(
    const argtype &args = argtype()) {
  return std::make_shared<TrialRemoveMultiple>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_REMOVE_MULTIPLE_H_
