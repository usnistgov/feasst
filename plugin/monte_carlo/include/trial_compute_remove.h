
#ifndef FEASST_MONTE_CARLO_TRIAL_COMPUTE_REMOVE_H_
#define FEASST_MONTE_CARLO_TRIAL_COMPUTE_REMOVE_H_

#include <vector>
#include "utils/include/arguments.h"
#include "system/include/system.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/trial_select.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute.h"

namespace feasst {

class TrialComputeRemove : public TrialCompute {
 public:
  TrialComputeRemove();

  void perturb_and_acceptance(
      Criteria * criteria,
      System * system,
      Acceptance * acceptance,
      std::vector<TrialStage*> * stages,
      Random * random) override {
    DEBUG("TrialComputeRemove");
    compute_rosenbluth(1, criteria, system, acceptance, stages, random);
    acceptance->set_energy_new(criteria->current_energy() - acceptance->energy_old());
    acceptance->add_to_macrostate_shift(-1);
    { // Metropolis
      const Configuration& config = system->configuration();
      const double volume = config.domain().volume();
      const TrialSelect * select = (*stages)[0]->trial_select();
      const int particle_index = select->mobile().particle_index(0);
      const int particle_type = config.select_particle(particle_index).type();
      DEBUG("volume " << volume << " selprob " << select->probability() << " betamu " << criteria->beta_mu(particle_type));
      acceptance->add_to_ln_metropolis_prob(
        - log(volume*select->probability())
        - criteria->beta_mu(particle_type)
      );
      DEBUG("lnmet " << acceptance->ln_metropolis_prob());
    }
  }
  std::shared_ptr<TrialCompute> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialComputeRemove(std::istream& istr);
  virtual ~TrialComputeRemove() {}

 protected:
  void serialize_trial_compute_remove_(std::ostream& ostr) const;
};

inline std::shared_ptr<TrialComputeRemove> MakeTrialComputeRemove() {
  return std::make_shared<TrialComputeRemove>();
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_COMPUTE_REMOVE_H_
