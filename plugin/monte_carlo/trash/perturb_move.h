#ifndef FEASST_MONTE_CARLO_PERTURB_MOVE_H_
#define FEASST_MONTE_CARLO_PERTURB_MOVE_H_

#include "monte_carlo/include/perturb.h"

namespace feasst {

class PerturbSelectMove : public PerturbOptRevert {
 public:
  void revert() override {
    DEBUG("revert possible " << revert_possible());
    if (revert_possible()) {
      Configuration* config = system()->get_configuration();
      DEBUG("reverting positions: " << selection().str());
      config->update_positions(selection(), false); // don't wrap if reverting
      system()->revert();
    }
  }

  Configuration * get_config_before_move(System * system) {
    store_old(system);
    // ASSERT(selection().num_sites() > 0, "no selection");
    return system->get_configuration();
  }

  void after_move() {
    set_revert_possible();
    set_selection_state("move");
  }

  ~PerturbSelectMove() {}
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_PERTURB_MOVE_H_
