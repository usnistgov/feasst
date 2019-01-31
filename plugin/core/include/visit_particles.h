
#ifndef FEASST_CORE_VISIT_PARTICLES_H_
#define FEASST_CORE_VISIT_PARTICLES_H_

#include "core/include/particles.h"
#include "core/include/select.h"

namespace feasst {

class LoopOneBody;

class VisitParticles {
 public:
  void loop(const Particles& particles,
            LoopOneBody * loop,
            const Select& select);
  virtual ~VisitParticles() {}
};

class LoopOneBody {
 public:
  void visit(VisitParticles& visitor,
             const Particles& particles,
             const Select& select) {
    visitor.loop(particles, this, select);
  }
  virtual void work(const Site& site) const = 0;
  virtual ~LoopOneBody() {}
};

}  // namespace feasst

#endif  // FEASST_CORE_VISIT_PARTICLES_H_
