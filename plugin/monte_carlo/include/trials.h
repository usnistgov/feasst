#ifndef FEASST_MONTE_CARLO_TRIALS_H_
#define FEASST_MONTE_CARLO_TRIALS_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_factory.h"
#include "monte_carlo/include/trial_move.h"

namespace feasst {

/// Attempt a rigid translation of a random particle.
class TrialTranslate : public TrialMove {
 public:
  TrialTranslate(argtype args = argtype());
  TrialTranslate(argtype * args);
  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialTranslate>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialTranslate>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialTranslate(std::istream& istr);
  virtual ~TrialTranslate() {}
};

inline std::shared_ptr<TrialTranslate> MakeTrialTranslate(argtype args = argtype()) {
  return std::make_shared<TrialTranslate>(args); }

/// Attempt a rigid rotation of a random particle.
std::shared_ptr<Trial> MakeTrialRotate(argtype args = argtype());

/// Attempt to add a particle.
class TrialAdd : public Trial {
 public:
  TrialAdd(argtype args = argtype());
  TrialAdd(argtype * args);
  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialAdd>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialAdd>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialAdd(std::istream& istr);
  virtual ~TrialAdd() {}
};

inline std::shared_ptr<TrialAdd> MakeTrialAdd(argtype args = argtype()) {
  return std::make_shared<TrialAdd>(args); }

/// Attempt to remove a particle.
std::shared_ptr<Trial> MakeTrialRemove(argtype args = argtype());

/// Attempt TrialAdd or TrialRemove with equal probability.
std::shared_ptr<TrialFactory> MakeTrialTransfer(
  argtype args = argtype());

/// Attempt to change the volume.
std::shared_ptr<Trial> MakeTrialVolume(argtype args = argtype());

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIALS_H_
