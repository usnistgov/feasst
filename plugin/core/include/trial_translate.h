
#ifndef FEASST_CORE_TRIAL_TRANSLATE_H_
#define FEASST_CORE_TRIAL_TRANSLATE_H_

#include "core/include/trial.h"
#include "core/include/perturb_translate.h"
#include "core/include/random.h"
#include "core/include/utils_io.h"

namespace feasst {

/**
 */
class TrialTranslate : public Trial {
 public:
  TrialTranslate() { set_group_index(); }

  void attempt(Criteria* criteria, System * system) {
    perturb_.before_attempt();
    criteria->before_attempt(system);
    perturb_.select_random_particle(group_index(), system->configuration());
    if (perturb_.selection().is_empty()) {
      // no particles present
      accept_criteria_.force_rejection = 1;
    } else {
      const double pe_old = system->energy(perturb_.selection());
      DEBUG("pe_old " << pe_old);
      const Position trajectory = random_.position_in_cube(
        system->dimension(),
        max_move_
      );
      perturb_.translate_selected_particle(trajectory, system);
      const double pe_new = system->energy(perturb_.selection());
      DEBUG("pe_new " << pe_new);
      const double delta_energy = pe_new - pe_old;
      accept_criteria_.ln_metropolis_prob = -criteria->beta()*delta_energy;
      accept_criteria_.energy_new = criteria->running_energy() + delta_energy;
      accept_criteria_.force_rejection = 0;
      accept_criteria_.system = system;
      DEBUG("delta_energy " << delta_energy);
    }
    if (criteria->is_accepted(accept_criteria_)) {
      record_success();
    } else {
      perturb_.revert();
    }
  }

  int verbose = 0;

  virtual ~TrialTranslate() {}

  void set_max_move(const double max_move) { max_move_ = max_move; }

 private:
  PerturbTranslate perturb_;
  double max_move_ = 0.1;
  Random random_;
  AcceptanceCriteria accept_criteria_;
};

}  // namespace feasst

#endif  // FEASST_CORE_TRIAL_TRANSLATE_H_
