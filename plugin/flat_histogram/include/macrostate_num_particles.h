
#ifndef FEASST_FLAT_HISTOGRAM_MACROSTATE_NUM_PARTICLES_H_
#define FEASST_FLAT_HISTOGRAM_MACROSTATE_NUM_PARTICLES_H_

#include "flat_histogram/include/macrostate.h"

namespace feasst {

/**
  Defines the macrostate to be the total number of particles in the system.
 */
class MacrostateNumParticles : public Macrostate {
 public:
  MacrostateNumParticles(const Histogram& histogram) : Macrostate(histogram) {}

  double value(const System* system, const Criteria* criteria) {
    return static_cast<double>(system->configuration().num_particles());
  }

  virtual ~MacrostateNumParticles() {}
};

inline std::shared_ptr<MacrostateNumParticles> MakeMacrostateNumParticles(
    const Histogram& histogram) {
  return std::make_shared<MacrostateNumParticles>(histogram);
}

}  // namespace feasst

#endif  // FEASST_FLAT_HISTOGRAM_MACROSTATE_NUM_PARTICLES_H_
