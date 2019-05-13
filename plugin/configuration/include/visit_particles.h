
#ifndef FEASST_CONFIGURATION_VISIT_PARTICLES_H_
#define FEASST_CONFIGURATION_VISIT_PARTICLES_H_

#include "configuration/include/particle_factory.h"
#include "configuration/include/select.h"

namespace feasst {

class LoopOneBody;

class VisitParticles {
 public:
  void loop(const ParticleFactory& particles,
            LoopOneBody * loop,
            const Select& select);
  virtual ~VisitParticles() {}
};

class LoopOneBody {
 public:
  void visit(VisitParticles& visitor,
             const ParticleFactory& particles,
             const Select& select) {
    visitor.loop(particles, this, select);
  }
  virtual void work(const Site& site) const = 0;
  virtual ~LoopOneBody() {}
};

}  // namespace feasst

#endif  // FEASST_CONFIGURATION_VISIT_PARTICLES_H_
