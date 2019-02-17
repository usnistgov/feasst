#include <vector>
#include "core/include/visit_particles.h"

namespace feasst {

void VisitParticles::loop(const ParticleFactory& particles,
                          LoopOneBody * loop,
                          const Select& select) {
  for (int select_index = 0;
       select_index < select.num_particles();
       ++select_index) {
    const int part_index = select.particle_index(select_index);
    const Particle part = particles.particle(part_index);
    for (int site_index : select.site_indices(select_index)) {
      const Site& site = part.sites()[site_index];
      loop->work(site);
    }
  }
}

}  // namespace feasst
