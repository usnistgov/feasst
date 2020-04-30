
#ifndef FEASST_MONTE_CARLO_TRIAL_TRANSLATE_H_
#define FEASST_MONTE_CARLO_TRIAL_TRANSLATE_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_move.h"

namespace feasst {

/**
  Attempt a rigid translation of a random particle.
 */
class TrialTranslate : public TrialMove {
 public:
  /// These arguments are sent to both PerturbTranslate and TrialStage.
  explicit TrialTranslate(const argtype& args = argtype());
  std::shared_ptr<Trial> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialTranslate(std::istream& istr);
  virtual ~TrialTranslate() {}

 protected:
  void serialize_trial_translate_(std::ostream& ostr) const;
};

inline std::shared_ptr<TrialTranslate> MakeTrialTranslate(
    const argtype &args = argtype()) {
  return std::make_shared<TrialTranslate>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_TRANSLATE_H_
