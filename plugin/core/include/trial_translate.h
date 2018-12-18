
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
    perturb_.select_random_particle(group_index(), system->config());
    if (perturb_.selection().is_empty()) {
      // no particles present
      accept_criteria_.force_rejection = 1;
    } else {
      const double pe_old = system->energy(perturb_.selection());
      const Position trajectory = random_.position_in_cube(
        system->dimension(),
        max_move_
      );
      perturb_.translate_selected_particle(trajectory, system);
      const double pe_new = system->energy(perturb_.selection());
      const double delta_energy = pe_new - pe_old;
      accept_criteria_.ln_metropolis_prob = -criteria->beta()*delta_energy;
      accept_criteria_.energy_new = criteria->running_energy() + delta_energy;
      accept_criteria_.force_rejection = 0;
      accept_criteria_.system = system;
      accept_criteria_.pair = system->visitor();
    }
    if (criteria->is_accepted(accept_criteria_)) {
      //cout << "accept " << delta_energy << endl;
    } else {
      //cout << "reject " << delta_energy << endl;
      perturb_.revert();
    }
    if (verbose==1) std::cout << "position " << str(system->particle(0).position().coord()) << std::endl;
    if (verbose==1) std::cout << "position " << str(system->particle(0).site(0).position().coord()) << std::endl;
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
