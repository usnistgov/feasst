#ifndef FEASST_MONTE_CARLO_TRIAL_TRANSLATE_H_
#define FEASST_MONTE_CARLO_TRIAL_TRANSLATE_H_

#include <map>
#include <string>
#include <memory>
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_move.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/// Attempt a rigid translation of a random particle.
class TrialTranslate : public TrialMove {
 public:
  //@{
  /** @name Arguments
    - Trial arguments.
    - TrialStage arguments.
    - TrialSelectParticle arguments.
    - PerturbTranslate arguments.
    - Tunable arguments.
   */
  explicit TrialTranslate(argtype args = argtype());
  explicit TrialTranslate(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialTranslate>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialTranslate>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialTranslate(std::istream& istr);
  virtual ~TrialTranslate() {}
  //@}
};

inline std::shared_ptr<TrialTranslate> MakeTrialTranslate(
    argtype args = argtype()) {
  return std::make_shared<TrialTranslate>(args); }

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_TRANSLATE_H_
