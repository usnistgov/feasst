
#ifndef FEASST_MONTE_CARLO_PERTURB_TRANSLATE_H_
#define FEASST_MONTE_CARLO_PERTURB_TRANSLATE_H_

#include "monte_carlo/include/perturb_move.h"

namespace feasst {

class PerturbTranslate : public PerturbSelectMove {
 public:
  void perturb(System * system) override {
    const Position trajectory = random_.position_in_cube(
      system->dimension(),
      tunable().value()
    );
    DEBUG("max move " << tunable().value());
    ASSERT(tunable().value() > NEAR_ZERO, "tunable is too small");
    translate_selection(trajectory, system);
  }

  void translate_selection(const Position &trajectory,
    System * system) {
    Configuration * config = get_config_before_move(system);
    displaced_ = selection();
    displace(trajectory, &displaced_);
    config->update_positions(displaced_);
    after_move();
  }

  void displace(const Position& displacement, SelectList * displaced) {
    for (int select_index = 0;
         select_index < displaced->num_particles();
         ++select_index) {
      Position displaced_part(displaced->particle_positions()[select_index]);
      displaced_part.add(displacement);
      displaced->set_particle_position(select_index, displaced_part);
      for (int site = 0;
           site < static_cast<int>(displaced->site_indices(select_index).size());
           ++site) {
        Position displaced_site(displaced->site_positions()[select_index][site]);
        displaced_site.add(displacement);
        displaced->set_site_position(select_index, site, displaced_site);
      }
    }
  }

  ~PerturbTranslate() {}

 private:
  Random random_;
  SelectList displaced_;
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_PERTURB_TRANSLATE_H_
