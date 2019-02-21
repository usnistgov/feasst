
#ifndef FEASST_CORE_PERTURB_TRANSLATE_H_
#define FEASST_CORE_PERTURB_TRANSLATE_H_

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
    return system->get_configuration();
  }

  void after_move() {
    set_revert_possible();
    set_selection_state("move");
  }

  ~PerturbSelectMove() {}
};

class PerturbTranslate : public PerturbSelectMove {
 public:
  void translate_selection(const Position &trajectory,
    System * system) {
    Configuration * config = get_config_before_move(system);
    config->displace_particles(selection(), trajectory);
    after_move();
  }
  ~PerturbTranslate() {}
};

class PerturbRotate : public PerturbSelectMove {
 public:
  PerturbRotate() { set_recenter(); }

  void rotate_selection(const Position& pivot,
      const RotationMatrix& rotation,
      System * system) {
    Configuration * config = get_config_before_move(system);
    SelectList rotated = selection();
    // for each site/particle in selection
    for (int select_index = 0;
         select_index < rotated.num_particles();
         ++select_index) {
      // rotate site positions
      for (int site = 0;
           site < static_cast<int>(rotated.site_indices(select_index).size());
           ++site) {
        Position position = rotated.site_positions()[select_index][site];
        rotate_(pivot, rotation, &position);
        rotated.set_site_position(select_index, site, position);
      }

      // rotate or recenter particle positions
      if (recenter_ == 0) {
        Position position = rotated.particle_positions()[select_index];
        rotate_(pivot, rotation, &position);
        rotated.set_particle_position(select_index, position);
      } else {
        const int part_index = rotated.particle_index(select_index);
        const Position center = config->select_particle(part_index).average_site_position();
        rotated.set_particle_position(select_index, center);
      }
    }
    config->update_positions(rotated);
    after_move();
  }
  ~PerturbRotate() {}

  /// By default, do not recenter the particle position based on sites.
  void set_recenter(const int recenter = 0) { recenter_ = recenter; }

 private:
  int recenter_ = 0;

  void rotate_(const Position& pivot, const RotationMatrix& rotation, Position * pos) const {
    pos->subtract(pivot);
    *pos = rotation.multiply(*pos);
    pos->add(pivot);
  }
};

}  // namespace feasst

#endif  // FEASST_CORE_PERTURB_TRANSLATE_H_
