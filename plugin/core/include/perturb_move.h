#ifndef FEASST_CORE_PERTURB_MOVE_H_
#define FEASST_CORE_PERTURB_MOVE_H_

#include "core/include/perturb.h"

namespace feasst {

class PerturbSelectMove : public PerturbOptRevert {
 public:
  void revert() override {
    if (revert_possible()) {
      Configuration* config = system()->get_configuration();
      config->update_positions(selection());
      system()->revert();
    }
  }

  Configuration * get_config_before_move(System * system) {
    store_old(system);
    ASSERT(selection().num_sites() > 0, "no selection");
    return system->get_configuration();
  }

  void after_move() {
    set_revert_possible();
    set_selection_state("move");
  }

  ~PerturbSelectMove() {}
};

}  // namespace feasst

#endif  // FEASST_CORE_PERTURB_MOVE_H_
