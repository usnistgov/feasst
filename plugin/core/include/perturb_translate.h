
#ifndef FEASST_CORE_PERTURB_TRANSLATE_H_
#define FEASST_CORE_PERTURB_TRANSLATE_H_

#include "core/include/perturb.h"

namespace feasst {

class PerturbTranslate : public PerturbOptRevert {
 public:
  const Select& selection() const override { return selection_; }

  void select_random_particle(const int group_index, const Configuration& config) {
    selection_.random_particle(config, group_index);
    selection_.set_trial_state("old");
  }

  void translate_selected_particle(const Position &trajectory,
    System * system) {
    store_old(system);
    Configuration* config = system->get_configuration();
    config->displace_particles(selection_, trajectory);
    set_revert_possible();
    selection_.set_trial_state("move");
  }

  void revert() override {
    if (revert_possible()) {
      Configuration* config = system()->get_configuration();
      config->update_positions(selection_);
      system()->revert();
    }
  }

  ~PerturbTranslate() {}

 private:
  SelectList selection_;
};

}  // namespace feasst

#endif  // FEASST_CORE_PERTURB_TRANSLATE_H_
