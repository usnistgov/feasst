
#ifndef FEASST_CORE_VISIT_CONFIGURATION_H_
#define FEASST_CORE_VISIT_CONFIGURATION_H_

#include "core/include/visit_particles.h"
#include "core/include/configuration.h"

namespace feasst {

class VisitConfiguration : public VisitParticles {
 public:
  virtual void loop(const Configuration& config,
                    LoopOneBody * loop_one_body,
                    const int group_index);
  virtual ~VisitConfiguration() {}
};

}  // namespace feasst

#endif  // FEASST_CORE_VISIT_CONFIGURATION_H_
