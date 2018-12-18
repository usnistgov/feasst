#include <vector>
#include "core/include/visit_configuration.h"

namespace feasst {

// HWH depreciate
void VisitConfiguration::loop(const Configuration& config,
                              LoopOneBody * loop_one_body,
                              const int group_index) {
  VisitParticles::loop(config.particles(), loop_one_body, config.group_selects()[group_index]);
}

}  // namespace feasst
