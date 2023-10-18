
#ifndef FEASST_CHARGE_TRIAL_REMOVE_MULTIPLE_H_
#define FEASST_CHARGE_TRIAL_REMOVE_MULTIPLE_H_

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
  //@{
  /** @name Arguments
   */

  /**
    args:
    - particle_type[i]: the i-th type of particle to add.
      The "[i]" is to be substituted for an integer 0, 1, 2, ...
    - TrialStage arguments.
   */
  explicit TrialRemoveMultiple(argtype args = argtype());
  explicit TrialRemoveMultiple(argtype * args);
  //@}
  /** @name Public Functions
   */
  //@{
  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialRemoveMultiple>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialRemoveMultiple>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialRemoveMultiple(std::istream& istr);
  virtual ~TrialRemoveMultiple() {}
  //@}
};

inline std::shared_ptr<TrialRemoveMultiple> MakeTrialRemoveMultiple(argtype args = argtype()) {
  return std::make_shared<TrialRemoveMultiple>(args); }

}  // namespace feasst

#endif  // FEASST_CHARGE_TRIAL_REMOVE_MULTIPLE_H_
