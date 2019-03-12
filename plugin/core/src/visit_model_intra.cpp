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
  TRACE("intra particle energy_of_selection");
  ASSERT(group_index == 0, "need to implement site1 loop filtering particles by group");
  zero_energy();
  const Domain& domain = config->domain();
  Position relative;
  relative.set_vector(domain.side_length().coord());
  for (int part1_index : selection.particle_indices()) {
    TRACE("particle: " << part1_index);
    const Particle& part1 = config->select_particle(part1_index);
    // the first site loop is over all sites in part1 and group_index
    // the second is all sites in selection
    // but this creates non-trivial double counting/optimization issues
    // for ease of implementation, we simply loop over all unique pairs of sites

    // here we use excluded to account for chain regrowth, etc.
    // exclude the particles which haven't been grown yet.
    Select sites;
    sites.add_particle(part1, part1_index);
    if (selection.excluded()) {
      sites.remove(*(selection.excluded()));
      TRACE("excluded " << selection.excluded()->str());
    }
    TRACE("sites: " << sites.str());
    const std::vector<int>& site_indices = sites.site_indices(0);
    for (int sel1_index = 0;
         sel1_index < static_cast<int>(site_indices.size()) - 1;
         ++sel1_index) {
      TRACE("sel1_index " << sel1_index << " size " << site_indices.size());
      const int site1_index = site_indices[sel1_index];
      for (int sel2_index = sel1_index + 1;
           sel2_index < static_cast<int>(site_indices.size());
           ++sel2_index) {
        const int site2_index = site_indices[sel2_index];
        if (std::abs(site1_index - site2_index) > intra_cut_) {
          TRACE("sites: " << site1_index << " " << site2_index);
          inner()->compute(part1_index, site1_index, part1_index, site2_index,
                           config, model_params, model, &relative);
        }
      }
    }
  }
  set_energy(inner()->energy());
}

}  // namespace feasst
