
#ifndef FEASST_MONTE_CARLO_TRIAL_ROTATE_H_
#define FEASST_MONTE_CARLO_TRIAL_ROTATE_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_move.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/perturb_rotate.h"

namespace feasst {

/// Attempt a rigid rotation of a random particle.
class TrialRotate : public TrialMove {
 public:
  TrialRotate(
    /// These arguments are sent to both PerturbRotate and TrialStage.
    const argtype& args = argtype())
    : TrialMove(
      std::make_shared<TrialSelectParticle>(),
      std::make_shared<PerturbRotate>(args),
      args
    ) {
    class_name_ = "TrialRotate";
  }
  std::shared_ptr<Trial> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialRotate(std::istream& istr);
  virtual ~TrialRotate() {}

 protected:
  void serialize_trial_rotate_(std::ostream& ostr) const;
};

inline std::shared_ptr<TrialRotate> MakeTrialRotate(
    const argtype &args = argtype()) {
  return std::make_shared<TrialRotate>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_ROTATE_H_
