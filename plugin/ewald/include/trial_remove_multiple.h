
#ifndef FEASST_MONTE_CARLO_TRIAL_REMOVE_MULTIPLE_H_
#define FEASST_MONTE_CARLO_TRIAL_REMOVE_MULTIPLE_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"

namespace feasst {

/**
  Attempt to remove multiple particles.
  Typically requires the use of a reference index.
 */
class TrialRemoveMultiple : public Trial {
 public:
  /**
    args:
    - particle_type[i]: the i-th type of particle to add.
      The "[i]" is to be substituted for an integer 0, 1, 2, ...
   */
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
