#ifndef FEASST_MONTE_CARLO_TRIAL_ROTATE_H_
#define FEASST_MONTE_CARLO_TRIAL_ROTATE_H_

#include <map>
#include <string>
#include <memory>
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_move.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
  Attempt a rigid rotation of a random particle.
  In most cases, TrialParticlePivot is more efficient.
 */
class TrialRotate : public TrialMove {
 public:
  //@{
  /** @name Arguments
   */

  /**
    args:
    - TrialSelectParticle arguments.
    - TrialStage arguments.
    - Tunable arguments.
   */
  explicit TrialRotate(argtype args = argtype());
  explicit TrialRotate(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialRotate>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialRotate>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialRotate(std::istream& istr);
  virtual ~TrialRotate() {}
  //@}
};

inline std::shared_ptr<TrialRotate> MakeTrialRotate(argtype args = argtype()) {
  return std::make_shared<TrialRotate>(args); }

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_ROTATE_H_
