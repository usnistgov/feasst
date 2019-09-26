
#ifndef FEASST_MONTE_CARLO_TRIAL_REMOVE_H_
#define FEASST_MONTE_CARLO_TRIAL_REMOVE_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/trial_compute_remove.h"
#include "monte_carlo/include/perturb_remove.h"

namespace feasst {

/// Attempt to remove a particle.
class TrialRemove : public Trial {
 public:
  TrialRemove(const argtype& args = argtype()) : Trial(args) {
    argtype args2(args);
    args2.insert({"load_coordinates", "false"});
    add_stage(
      std::make_shared<TrialSelectParticle>(args2),
      std::make_shared<PerturbRemove>(),
      args
    );
    set(std::make_shared<TrialComputeRemove>());
    class_name_ = "TrialRemove";
  }
  std::shared_ptr<Trial> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialRemove(std::istream& istr);
  virtual ~TrialRemove() {}

 protected:
  void serialize_trial_remove_(std::ostream& ostr) const;
};

inline std::shared_ptr<TrialRemove> MakeTrialRemove(
    const argtype &args = argtype()) {
  return std::make_shared<TrialRemove>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_REMOVE_H_
