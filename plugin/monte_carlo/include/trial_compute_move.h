
#ifndef FEASST_MONTE_CARLO_TRIAL_COMPUTE_MOVE_H_
#define FEASST_MONTE_CARLO_TRIAL_COMPUTE_MOVE_H_

#include <vector>
#include "utils/include/arguments.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute.h"

namespace feasst {

class TrialComputeMove : public TrialCompute {
 public:
  TrialComputeMove();

  void perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) override {
    DEBUG("TrialComputeMove");
    compute_rosenbluth(1, criteria, system, acceptance, stages, random);
    for (TrialStage * stage : *stages) stage->mid_stage(system);
    compute_rosenbluth(0, criteria, system, acceptance, stages, random);
    DEBUG("old: " << criteria->current_energy() << " " << acceptance->energy_old());
    DEBUG("new: " << acceptance->energy_new());
    DEBUG("energy change: " << acceptance->energy_new() - acceptance->energy_old());
    const double delta_energy = acceptance->energy_new() - acceptance->energy_old();
    acceptance->set_energy_new(criteria->current_energy() + delta_energy);
  }
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialComputeMove(std::istream& istr);
  virtual ~TrialComputeMove() {}

 protected:
  void serialize_trial_compute_move_(std::ostream& ostr) const;
};

inline std::shared_ptr<TrialComputeMove> MakeTrialComputeMove() {
  return std::make_shared<TrialComputeMove>();
}
}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_COMPUTE_MOVE_H_
