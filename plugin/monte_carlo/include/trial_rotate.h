#ifndef FEASST_MONTE_CARLO_TRIAL_ROTATE_H_
#define FEASST_MONTE_CARLO_TRIAL_ROTATE_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_move.h"

namespace feasst {

/// Attempt a rigid rotation of a random particle.
class TrialRotate : public TrialMove {
 public:
  TrialRotate(argtype args = argtype());
  TrialRotate(argtype * args);
  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialRotate>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialRotate>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialRotate(std::istream& istr);
  virtual ~TrialRotate() {}
};

inline std::shared_ptr<TrialRotate> MakeTrialRotate(argtype args = argtype()) {
  return std::make_shared<TrialRotate>(args); }

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_ROTATE_H_
