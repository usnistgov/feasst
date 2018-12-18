
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
  }

  void add(const Particle &particle, System * system) {
    store_old(system);
    added_ = 1;
    Configuration* config = system->configuration(0);
    config->add_particle(particle.type());
    selection_.last_particle_added(config);
    config->replace_position(selection_, particle);
    set_revert_possible();
  }

  void remove_selected_particle(System * system) {
    store_old(system);
    added_ = 0;
    Configuration* config = system->configuration(0);
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
          system()->configuration(0)->remove_particles(selection_);
        } else {
          TRACE("reverting deletion");
          system()->configuration(0)->revive(selection_);
        }
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
