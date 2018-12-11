
#ifndef FEASST_CORE_PERTURB_TRANSLATE_H_
#define FEASST_CORE_PERTURB_TRANSLATE_H_

#include "core/include/perturb.h"

namespace feasst {

class PerturbTranslate : public Perturb {
 public:
  void translate_selected_particle(const Position &trajectory,
    System * system) {
    store_old(system);
    Configuration* config = system->configuration(0);
    if (optimization_ != 0) {
      particle_ = config->selected_particle();
      selection_ = config->selection();
    }
    config->displace_selected_particles(trajectory);
    set_revert_possible();
  }

  void revert() {
    if (revert_possible()) {
      if (optimization_ == 0) {
        Perturb::revert();
      } else {
        Configuration* config = system()->configuration(0);
        config->set_selection(selection_);
        config->replace_selected_particle_position(particle_);
      }
    }
  }

  ~PerturbTranslate() {}

 private:
  Selection selection_;
  Particle particle_;
};

}  // namespace feasst

#endif  // FEASST_CORE_PERTURB_TRANSLATE_H_
