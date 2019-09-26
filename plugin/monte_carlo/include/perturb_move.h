
#ifndef FEASST_MONTE_CARLO_PERTURB_MOVE_H_
#define FEASST_MONTE_CARLO_PERTURB_MOVE_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/perturb.h"

namespace feasst {

/**
  Only perturb the positions of the particles and/or sites.
 */
class PerturbMove : public Perturb {
 public:
  PerturbMove(const argtype& args = argtype()) : Perturb(args) {}

  /// Move the selection of the system.
  virtual void move(System * system, TrialSelect * select, Random * random) = 0;

  // The perturbation move is simplified such that the move of the selection of
  // the system is all that remains to be implemented.
  void perturb(
      System * system,
      TrialSelect * select,
      Random * random,
      const bool is_position_held = false
      ) override {
    if (is_position_held) {
      select->set_trial_state(0);
      return;
    }
    move(system, select, random);
    set_revert_possible(true, select);
    select->set_trial_state(1);
  }

  /// For perturbations that only move particles and/or sites, the revert step
  /// is the same for all. Simply put the original positions back.
  void revert(System * system) override {
    if (revert_possible()) {
      Configuration* config = system->get_configuration();
      config->update_positions(revert_select()->mobile_original(),
        // don't wrap if reverting
        false);
      system->revert();
    }
  }

  PerturbMove(std::istream& istr) : Perturb(istr) {}
  virtual ~PerturbMove() {}
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_PERTURB_MOVE_H_
