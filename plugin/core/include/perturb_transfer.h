
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

  const Select& selection() const override { return selection_; }

  void select_random_particle(const int type, Configuration * config) {
    ASSERT(optimization_ == 1, "error");
    const int load_coordinates = 0;  // don't load coordinates
    selection_.random_particle_of_type(type, config, load_coordinates);
    selection_.set_trial_state("old");
  }

  void add(const Particle &particle, System * system) {
    store_old(system);
    added_ = 1;
    Configuration* config = system->get_configuration();
    config->add_particle(particle.type());
    selection_.last_particle_added(config);
    selection_.set_trial_state("add");
    config->replace_position(selection_, particle);
    set_revert_possible();
  }

  void remove_selected_particle(System * system) {
    store_old(system);
    added_ = 0;
    Configuration* config = system->get_configuration();
    if (optimization_ != 0) {
      config->remove_particles(selection_);
    }
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
          TRACE("reverting addition");
          system()->get_configuration()->remove_particles(selection_);
        } else {
          TRACE("reverting deletion");
          system()->get_configuration()->revive(selection_);
        }
        system()->revert();
      }
      TRACE("done reverting");
    }
  }

  ~PerturbTransfer() {}

 private:
  int added_;
  SelectList selection_;
};

}  // namespace feasst

#endif  // FEASST_CORE_PERTURB_TRANSFER_H_
