
#ifndef FEASST_MONTE_CARLO_TRIAL_TRANSLATE_H_
#define FEASST_MONTE_CARLO_TRIAL_TRANSLATE_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_move.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/perturb_translate.h"

namespace feasst {

/// Attempt a rigid translation of a random particle.
class TrialTranslate : public TrialMove {
 public:
  TrialTranslate(
    /// These arguments are sent to both PerturbTranslate and TrialStage.
    const argtype& args = argtype())
    : TrialMove(
      std::make_shared<TrialSelectParticle>(),
      std::make_shared<PerturbTranslate>(args),
      args
    ) {
    class_name_ = "TrialTranslate";
  }
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
