
#ifndef FEASST_CONFINEMENT_TRIAL_ANYWHERE_NEW_ONLY_H_
#define FEASST_CONFINEMENT_TRIAL_ANYWHERE_NEW_ONLY_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_move.h"

namespace feasst {

/**
  Attempt to rigidly move anywhere in the box with any orientation.
  Do not compute the energy of the old configuration.
 */
class TrialAnywhereNewOnly : public TrialMove {
 public:
  /// These arguments are sent to both PerturbAnywhere and TrialStage.
  explicit TrialAnywhereNewOnly(const argtype& args = argtype());
  std::shared_ptr<Trial> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialAnywhereNewOnly(std::istream& istr);
  virtual ~TrialAnywhereNewOnly() {}

 protected:
  void serialize_trial_anywhere_new_only_(std::ostream& ostr) const;
};

inline std::shared_ptr<TrialAnywhereNewOnly> MakeTrialAnywhereNewOnly(
    const argtype &args = argtype()) {
  return std::make_shared<TrialAnywhereNewOnly>(args);
}

}  // namespace feasst

#endif  // FEASST_CONFINEMENT_TRIAL_ANYWHERE_NEW_ONLY_H_
