
#ifndef FEASST_MONTE_CARLO_TRIAL_ADD_H_
#define FEASST_MONTE_CARLO_TRIAL_ADD_H_

#include <memory>
#include "utils/include/arguments.h"
#include "system/include/system.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/trial_compute_add.h"
#include "monte_carlo/include/perturb_add.h"

namespace feasst {

/// Attempt to add a particle.
class TrialAdd : public Trial {
 public:
  TrialAdd(const argtype& args = argtype()) : Trial(args) {
    add_stage(
      std::make_shared<TrialSelectParticle>(args),
      std::make_shared<PerturbAdd>(args),
      args
    );
    set(std::make_shared<TrialComputeAdd>());
    class_name_ = "TrialAdd";
  }
  std::shared_ptr<Trial> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialAdd(std::istream& istr);
  virtual ~TrialAdd() {}

 protected:
  void serialize_trial_add_(std::ostream& ostr) const;
};

inline std::shared_ptr<TrialAdd> MakeTrialAdd(
    const argtype &args = argtype()) {
  return std::make_shared<TrialAdd>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_ADD_H_
