#include <vector>
#include "core/include/visit_model_intra.h"
#include "core/include/model_two_body.h"

namespace feasst {

void VisitModelIntra::compute(
    const ModelTwoBody& model,
    const ModelParams& model_params,
    const Select& selection,
    Configuration * config,
    const int group_index) {
  DEBUG("intra particle energy_of_selection");
  ASSERT(group_index == 0, "need to implement site1 loop filtering particles by group");
  set_energy(0);
  const Domain& domain = config->domain();
  Position relative;
  relative.set_vector(domain.side_length().coord());
  for (int part1_index : selection.particle_indices()) {
    TRACE("particle: " << part1_index);
    const Particle part1 = config->select_particle(part1_index);
    // the first site loop is over all sites in part1 and group_index
    // the second is all sites in selection
    // but this creates non-trivial double counting/optimization issues
    // for ease of implementation, we simply loop over all unique pairs of sites
    for (int site1_index = 0;
         site1_index < part1.num_sites() - 1;
         ++site1_index) {
      const Site site1 = part1.sites()[site1_index];
      for (int site2_index = site1_index + 1;
           site2_index < part1.num_sites();
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
