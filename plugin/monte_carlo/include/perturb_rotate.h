
#ifndef FEASST_MONTE_CARLO_PERTURB_ROTATE_H_
#define FEASST_MONTE_CARLO_PERTURB_ROTATE_H_

#include "monte_carlo/include/perturb_move.h"

namespace feasst {

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
  /// Recentering uses the old configuration, not the new one, so the particles
  /// will not be perfectly recentered.
  void set_recenter(const int recenter = 0) { recenter_ = recenter; }

 private:
  int recenter_ = 0;

  void rotate_(const Position& pivot, const RotationMatrix& rotation, Position * pos) const {
    TRACE("rotating pivot " << pivot.str() << " pos " << pos->str());
    TRACE("matrix " << rotation.str());
    pos->subtract(pivot);
    *pos = rotation.multiply(*pos);
    pos->add(pivot);
    TRACE("new pos " << pos->str());
  }
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_PERTURB_ROTATE_H_
