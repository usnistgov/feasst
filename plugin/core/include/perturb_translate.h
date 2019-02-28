
#ifndef FEASST_CORE_PERTURB_TRANSLATE_H_
#define FEASST_CORE_PERTURB_TRANSLATE_H_

#include "core/include/perturb_move.h"

namespace feasst {

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

}  // namespace feasst

#endif  // FEASST_CORE_PERTURB_TRANSLATE_H_
