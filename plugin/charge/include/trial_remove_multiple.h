
#ifndef FEASST_CHARGE_TRIAL_REMOVE_MULTIPLE_H_
#define FEASST_CHARGE_TRIAL_REMOVE_MULTIPLE_H_

#include <memory>
#include "monte_carlo/include/trial.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
  Attempt to remove multiple particles.
  Typically requires the use of a reference index.
  For a derivation of the acceptance criteria, see TrialAddMultiple that is
  the reverse of this trial.
 */
class TrialRemoveMultiple : public Trial {
 public:
  //@{
  /** @name Arguments
    - particle_types: comma-separated list of particle names to remove.
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
