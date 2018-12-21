
#ifndef FEASST_CORE_MACROSTATE_NUM_PARTICLES_H_
#define FEASST_CORE_MACROSTATE_NUM_PARTICLES_H_

#include "core/include/macrostate.h"

namespace feasst {

/**
  Defines the macrostate to be the total number of particles in the system.
 */
class MacrostateNumParticles : public Macrostate {
 public:
  double value(const System* system, const Criteria* criteria) {
    return static_cast<double>(system->configuration().num_particles());
  }

  virtual ~MacrostateNumParticles() {}
};

}  // namespace feasst

#endif  // FEASST_CORE_MACROSTATE_NUM_PARTICLES_H_
