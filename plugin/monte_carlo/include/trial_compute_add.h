
#ifndef FEASST_MONTE_CARLO_TRIAL_COMPUTE_ADD_H_
#define FEASST_MONTE_CARLO_TRIAL_COMPUTE_ADD_H_

#include <vector>
#include "utils/include/arguments.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_select.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute.h"

namespace feasst {

class TrialComputeAdd : public TrialCompute {
 public:
  TrialComputeAdd();

  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override {
    DEBUG("TrialComputeAdd");
    compute_rosenbluth(0, criteria, system, acceptance, stages, random);
    const TrialSelect * select = (*stages)[0]->trial_select();
    system->get_configuration()->revive(select->mobile());
    acceptance->set_energy_new(criteria->current_energy() + acceptance->energy_new());
    { // Metropolis
      const Configuration& config = system->configuration();
      const double volume = config.domain().volume();
      const int particle_index = select->mobile().particle_index(0);
      const int particle_type = config.select_particle(particle_index).type();
      DEBUG("volume " << volume << " selprob " << select->probability() << " betamu " << criteria->beta_mu(particle_type));
      acceptance->add_to_ln_metropolis_prob(
        log(volume*select->probability())
        + criteria->beta_mu(particle_type)
      );
    }
  }
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialComputeAdd(std::istream& istr);
  virtual ~TrialComputeAdd() {}

 protected:
  void serialize_trial_compute_add_(std::ostream& ostr) const;
};

inline std::shared_ptr<TrialComputeAdd> MakeTrialComputeAdd() {
  return std::make_shared<TrialComputeAdd>();
}
}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_COMPUTE_ADD_H_
