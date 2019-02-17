
#ifndef FEASST_CORE_PERTURB_TRANSFER_H_
#define FEASST_CORE_PERTURB_TRANSFER_H_

#include "core/include/perturb.h"

namespace feasst {

class PerturbAdd : public PerturbOptRevert {
 public:
  const Select& selection() const override { return selection_; }

  void add(const int particle_type, const Position& position, System * system) {
    store_old(system);
    Configuration* config = system->get_configuration();
    config->add_particle(particle_type);
    selection_.last_particle_added(config);
    selection_.set_trial_state("add");
    config->displace_particles(selection_, position);
    set_revert_possible();
  }

  void revert() override {
    if (revert_possible()) {
      TRACE("reverting addition");
      system()->get_configuration()->remove_particles(selection_);
      system()->revert();
      TRACE("done reverting");
    }
  }

  ~PerturbAdd() {}

 private:
  SelectList selection_;
};

class PerturbRemove : public PerturbOptRevert {
 public:
  const Select& selection() const override { return selection_; }

  void select_random_particle(const int type, Configuration * config) {
    const int load_coordinates = 0;  // don't load coordinates
    selection_.random_particle_of_type(type, config, load_coordinates);
    selection_.set_trial_state("old");
  }

  void remove_selected_particle(System * system) {
    store_old(system);
    Configuration* config = system->get_configuration();
    config->remove_particles(selection_);
    set_revert_possible();
  }

  void revert() override {
    if (revert_possible()) {
      TRACE("reverting deletion");
      system()->get_configuration()->revive(selection_);
      system()->revert();
      TRACE("done reverting");
    }
  }

  ~PerturbRemove() {}

 private:
  SelectList selection_;
};

}  // namespace feasst

#endif  // FEASST_CORE_PERTURB_TRANSFER_H_
