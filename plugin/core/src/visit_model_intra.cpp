#include <vector>
#include "core/include/visit_model_intra.h"
#include "core/include/model_two_body.h"

namespace feasst {

void VisitModelIntra::compute(const Configuration& config,
                              const ModelTwoBody& model,
                              const Select& selection,
                              const int group_index) {
  DEBUG("intra particle energy_of_selection");
  double energy = 0;
  Position difference;
  const Domain& domain = config.domain();
  const ModelParams& model_params = config.unique_types().model_params();
  for (int part1_index : selection.particle_indices()) {
    TRACE("particle: " << part1_index);
    const Particle part1 = config.particle(part1_index);
    for (int site1_index : selection.site_indices(part1_index)) {
      const Site site1 = part1.sites()[site1_index];
      for (int site2_index : selection.site_indices(part1_index)) {
        if (std::abs(site1_index - site2_index) > intra_cut_) {
          TRACE("sites: " << site1_index << " " << site2_index);
          const Site site2 = part1.sites()[site2_index];
          difference = site1.position();
          difference.subtract(site2.position());
          domain.wrap(&difference);
          energy += model.evaluate(difference, site1, site2, model_params);
          TRACE("energy: " << energy);
        }
      }
    }
  }
  set_energy(energy);
}

}  // namespace feasst
