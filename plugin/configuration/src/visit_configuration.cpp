#include <vector>
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/configuration.h"
#include "configuration/include/visit_configuration.h"

namespace feasst {

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
      const Site& site = part.sites()[site_index];
      if (site.is_physical()) {
        data_.site_index = site_index;
        loop->work(site, config, data_);
      }
    }
  }
}

void VisitConfiguration::loop(const Configuration& config,
          LoopConfigOneBody * loop_config_one_body,
          const int group_index) {
  loop(config, loop_config_one_body, config.group_select(group_index));
}

}  // namespace feasst
