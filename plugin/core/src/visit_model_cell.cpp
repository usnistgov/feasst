#include <vector>
#include <sstream>
#include "core/include/visit_model_cell.h"
#include "core/include/model_two_body.h"
#include "core/include/model_one_body.h"
#include "core/include/utils_math.h"

namespace feasst {

void VisitModelCell::compute(
    const ModelTwoBody& model,
    const ModelParams& model_params,
    Configuration * config,
    const int cell_index) {
  zero_energy();
  const Domain& domain = config->domain();
  ASSERT(cell_index < static_cast<int>(domain.cells().size()), "index error");
  const Cells& cells = domain.cells()[cell_index];
  Position relative;
  relative.set_vector(domain.side_length().coord());

  /** Loop index nomenclature
    ends in 1 or 2 to represent the pair
    cell -> id of cell
    select -> selection inside each cell
    select_index -> index of selection
    part -> particle
    part_index -> index of particle
    site -> site
   */

  // loop through neighboring cells where cell1 < cell2 only
  for (int cell1 = 0; cell1 < cells.num_total(); ++cell1) {
    const Select& select1 = cells.particles()[cell1];
    for (int cell2 : cells.neighbor()[cell1]) {
      if (cell1 < cell2) {
        const Select& select2 = cells.particles()[cell2];
        for (int select1_index = 0;
             select1_index < select1.num_particles();
             ++select1_index) {
          const int part1_index = select1.particle_index(select1_index);
          for (int select2_index = 0;
               select2_index < select2.num_particles();
               ++select2_index) {
            const int part2_index = select2.particle_index(select2_index);
            if (part1_index != part2_index) {
              for (int site1_index : select1.site_indices(select1_index)) {
                for (int site2_index : select2.site_indices(select2_index)) {
                  inner()->compute(part1_index, site1_index, part2_index, site2_index,
                                   config, model_params, model, &relative);
                }
              }
            }
          }
        }
      }
    }
  }

  // loop through the same cell only
  for (int cell1 = 0; cell1 < cells.num_total(); ++cell1) {
    const Select& select = cells.particles()[cell1];
    for (int select1_index = 0;
         select1_index < select.num_particles() - 1;
         ++select1_index) {
      const int part1_index = select.particle_index(select1_index);
      for (int select2_index = select1_index + 1;
           select2_index < select.num_particles();
           ++select2_index) {
        const int part2_index = select.particle_index(select2_index);
        if (part1_index != part2_index) {
          for (int site1_index : select.site_indices(select1_index)) {
            for (int site2_index : select.site_indices(select2_index)) {
              inner()->compute(part1_index, site1_index, part2_index, site2_index,
                               config, model_params, model, &relative);
            }
          }
        }
      }
    }
  }
  set_energy(inner()->energy());
}

void VisitModelCell::compute(
    const ModelTwoBody& model,
    const ModelParams& model_params,
    const Select& selection,
    Configuration * config,
    const int cell_index) {
  DEBUG("visiting model");
  zero_energy();
  const Domain& domain = config->domain();
  ASSERT(cell_index < static_cast<int>(domain.cells().size()), "were cells not initialized?");
  const Cells& cells = domain.cells()[cell_index];
  Position relative;
  relative.set_vector(domain.side_length().coord());
  std::stringstream ss;
  ss << "cell" << cell_index;
  const std::string cell_label = ss.str();
  for (int select1_index = 0;
       select1_index < selection.num_particles();
       ++select1_index) {
    const int part1_index = selection.particle_index(select1_index);
    const Particle& part1 = config->select_particle(part1_index);
    for (int site1_index : selection.site_indices(select1_index)) {
      const Site& site1 = part1.site(site1_index);
      const int cell1_index = feasst::round(site1.property(cell_label));
      for (int cell2_index : cells.neighbor()[cell1_index]) {
        const Select& cell2_parts = cells.particles()[cell2_index];
        for (int select2_index = 0;
             select2_index < cell2_parts.num_particles();
             ++select2_index) {
          const int part2_index = cell2_parts.particle_index(select2_index);
          if (part1_index != part2_index) {
            TRACE("indices " <<
                  feasst_str(cell2_parts.site_indices(select2_index)));
            for (int site2_index : cell2_parts.site_indices(select2_index)) {
              TRACE("index: " << part1_index << " " << part2_index << " " <<
                   site1_index << " " << site2_index);
              inner()->compute(part1_index, site1_index, part2_index, site2_index,
                               config, model_params, model, &relative);
            }
          }
        }
      }
    }
  }
  set_energy(inner()->energy());
}

}  // namespace feasst
