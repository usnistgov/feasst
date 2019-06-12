
#ifndef FEASST_MONTE_CARLO_PERTURB_TRANSFER_H_
#define FEASST_MONTE_CARLO_PERTURB_TRANSFER_H_

#include "monte_carlo/include/perturb.h"
#include "monte_carlo/include/perturb_translate.h"
#include "monte_carlo/include/perturb_rotate.h"

namespace feasst {

//class PerturbAdd : public PerturbTranslate {
class PerturbAdd : public PerturbOptRevert {
 public:
  void add(const int particle_type, const Position& position, System * system) {
    store_old(system);
    Configuration* config = system->get_configuration();
    config->add_particle_of_type(particle_type);
    select_last_particle_added(config);
    set_selection_state("add");
    SelectList new_(selection());
    const Position& pivot = new_.particle_positions()[0];
    rotate_.rotate(
      pivot,
      random_.rotation(pivot.dimension()),
      &new_);
    translate_.displace(position, &new_);
    config->update_positions(new_);
    set_revert_possible();
  }

  void revert() override {
    if (revert_possible()) {
      TRACE("reverting addition");
      system()->get_configuration()->remove_particles(selection());
      system()->revert();
      TRACE("done reverting");
    }
  }

  ~PerturbAdd() {}

 private:
  PerturbTranslate translate_;
  PerturbRotate rotate_;
  Random random_;
};

class PerturbRemove : public PerturbOptRevert {
 public:
  void select_random_particle(const int type, Configuration * config) {
    select_random_particle_of_type(type, config);
    set_selection_state("old");
  }

  void remove_selected_particle(System * system) {
    store_old(system);
    Configuration* config = system->get_configuration();
    config->remove_particles(selection());
    set_revert_possible();
  }

  void revert() override {
    if (revert_possible()) {
      TRACE("reverting deletion");
      system()->get_configuration()->revive(selection());
      system()->revert();
      TRACE("done reverting");
    }
  }

  ~PerturbRemove() {}
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_PERTURB_TRANSFER_H_
