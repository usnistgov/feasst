#include <vector>
#include "configuration/include/visit_configuration.h"

namespace feasst {

// HWH depreciate
void VisitConfiguration::loop(const Configuration& config,
                              LoopConfigOneBody * loop,
                              const Select& select) {
  for (int select_index = 0;
       select_index < select.num_particles();
       ++select_index) {
    const int part_index = select.particle_index(select_index);
    data_.particle_index = part_index;
    const Particle& part = config.select_particle(part_index);
    data_.particle_type = part.type();
    for (int site_index : select.site_indices(select_index)) {
      data_.site_index = site_index;
      const Site& site = part.sites()[site_index];
      loop->work(site, config, data_);
    }
  }
}

}  // namespace feasst
