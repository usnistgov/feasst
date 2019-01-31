#include <vector>
#include "core/include/visit_model_intra.h"
#include "core/include/model_two_body.h"

namespace feasst {

void VisitModelIntra::compute(
    const ModelTwoBody& model,
    const Select& selection,
    Configuration * config,
    const int group_index) {
  DEBUG("intra particle energy_of_selection");
  set_energy(0);
  const Domain& domain = config->domain();
  Position relative;
  relative.set_vector(domain.side_length().coord());
  const ModelParams& model_params = config->unique_types().model_params();
  for (int part1_index : selection.particle_indices()) {
    TRACE("particle: " << part1_index);
    const Particle part1 = config->particle(part1_index);
    const std::vector<int>& site_indices = selection.site_indices(part1_index);
    for (int site1_index = 0;
         site1_index < static_cast<int>(site_indices.size()) - 1;
         ++site1_index) {
      const Site site1 = part1.sites()[site1_index];
      for (int site2_index = site1_index + 1;
           site2_index < static_cast<int>(site_indices.size());
           ++site2_index) {
        if (std::abs(site1_index - site2_index) > intra_cut_) {
          TRACE("sites: " << site1_index << " " << site2_index);
          const Site site2 = part1.sites()[site2_index];
          inner_(site1, site2, domain, model_params, model, &relative);
        }
      }
    }
  }
}

}  // namespace feasst
