
#ifndef FEASST_SYSTEM_VISIT_CONFIGURATION_H_
#define FEASST_SYSTEM_VISIT_CONFIGURATION_H_

#include "configuration/include/visit_particles.h"
#include "configuration/include/configuration.h"

namespace feasst {

class VisitConfiguration : public VisitParticles {
 public:
  void loop(const Configuration& config,
            LoopOneBody * loop_one_body,
            const int group_index);
  virtual ~VisitConfiguration() {}
};

}  // namespace feasst

#endif  // FEASST_SYSTEM_VISIT_CONFIGURATION_H_
