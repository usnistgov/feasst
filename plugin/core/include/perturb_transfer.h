
#ifndef FEASST_CORE_PERTURB_TRANSFER_H_
#define FEASST_CORE_PERTURB_TRANSFER_H_

#include "core/include/perturb.h"

namespace feasst {

class PerturbTransfer : public Perturb {
 public:
  void before_attempt() override {
    added_ = -1;
    Perturb::before_attempt();
  }

  void add(const Particle &particle, System * system) {
    store_old(system);
    added_ = 1;
    Configuration* config = system->configuration(0);
    config->add_particle(particle.type());
    config->select_last_particle();
    if (optimization_ != 0) {
      selection_ = config->selection();
    }
    config->replace_position_of_last(particle);
    set_revert_possible();
  }

  void remove_selected_particle(System * system) {
    store_old(system);
    added_ = 0;
    Configuration* config = system->configuration(0);
    if (optimization_ != 0) {
      particle_ = config->selected_particle();
    }
    config->remove_selected_particles();
    set_revert_possible();
  }

  void revert() override {
    if (revert_possible()) {
      if (optimization_ == 0) {
        Perturb::revert();
      } else {
        ASSERT(added_ == 1 || added_ == 0,
          "unrecognized added(" << added_ << ")");
        if (added_ == 1) {
          system()->configuration(0)->set_selection(selection_);
          remove_selected_particle(system());
        } else {
          add(particle_, system());
        }
      }
    }
  }

  ~PerturbTransfer() {}

 private:
  int added_;
  Particle particle_;
  Selection selection_;
};

}  // namespace feasst

#endif  // FEASST_CORE_PERTURB_TRANSFER_H_
