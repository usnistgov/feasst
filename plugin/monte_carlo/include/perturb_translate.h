
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
    config->displace_particles(selection(), trajectory);
    after_move();
  }
  ~PerturbTranslate() {}

 private:
  Random random_;
};

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_PERTURB_TRANSLATE_H_
