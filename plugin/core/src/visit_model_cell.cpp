#include <vector>
#include <sstream>
#include "core/include/visit_model_cell.h"
#include "core/include/model_two_body.h"
#include "core/include/model_one_body.h"
#include "core/include/utils_math.h"

namespace feasst {

void VisitModelCell::compute(const Configuration& config,
    const ModelTwoBody& model,
    const int cell_index) {
  double energy = 0;
  const Domain& domain = config.domain();
  ASSERT(cell_index < static_cast<int>(domain.cells().size()), "index error");
  const Cells& cells = domain.cells()[cell_index];
  Position relative;
  double r2 = 0.;
  relative.set_vector(domain.side_length().coord());
  const ModelParams& model_params = config.unique_types().model_params();

  /** Loop index nomenclature
    ends in 1 or 2 to represent the pair
    cell -> id of cell
    select -> selection inside each cell
    select_index -> index of selection
    part -> particle
    part_index -> index of particle
    site -> site
   */

  for (int cell1 = 0; cell1 < cells.num_total(); ++cell1) {
    const Select& select1 = cells.particles()[cell1];
    for (int cell2 : cells.neighbor()[cell1]) {
      const Select& select2 = cells.particles()[cell2];
      for (int select1_index = 0;
           select1_index < select1.num_particles();
           ++select1_index) {
        const int part1_index = select1.particle_index(select1_index);
        const Particle part1 = config.select_particle(part1_index);
        for (int select2_index = 0;
             select2_index < select2.num_particles();
             ++select2_index) {
          const int part2_index = select2.particle_index(select2_index);
          if (part1_index < part2_index) {
            const Particle part2 = config.select_particle(part2_index);
            for (int site1_index : select1.site_indices(select1_index)) {
              const Site& site1 = part1.site(site1_index);
              for (int site2_index : select2.site_indices(select2_index)) {
                const Site& site2 = part2.site(site2_index);
                domain.wrap_opt(site1.position(), site2.position(), &relative, &r2);
                energy += model.evaluate(relative, site1, site2, model_params);
              }
            }
          }
        }
      }
    }
  }
  set_energy(energy);
}

void VisitModelCell::compute(const Configuration& config,
                         const ModelTwoBody& model,
                         const Select& selection,
                         const int cell_index) {
  double energy = 0;
  const Domain& domain = config.domain();
  ASSERT(cell_index < static_cast<int>(domain.cells().size()), "index error");
  const Cells& cells = domain.cells()[cell_index];
  Position relative;
  double r2 = 0.;
  relative.set_vector(domain.side_length().coord());
  std::stringstream ss;
  ss << "cell" << cell_index;
  const std::string cell_label = ss.str();
  const ModelParams& model_params = config.unique_types().model_params();
  for (int select_index = 0; select_index < selection.num_particles(); ++select_index) {
    const int part1_index = selection.particle_index(select_index);
    const Particle part1 = config.select_particle(part1_index);
    for (int site1_index : selection.site_indices(select_index)) {
      const Site& site1 = part1.sites()[site1_index];
      const int cell1_index = feasst::round(site1.property(cell_label));
      for (int cell2_index : cells.neighbor()[cell1_index]) {
        const Select& cell2_parts = cells.particles()[cell2_index];
        for (int part2_index : cell2_parts.particle_indices()) {
          if (part1_index != part2_index) {
            const Particle part2 = config.select_particle(part2_index);
            for (int site2_index : cell2_parts.site_indices(select_index)) {
              const Site& site2 = part2.sites()[site2_index];
              domain.wrap_opt(site1.position(), site2.position(), &relative, &r2);
              energy += model.evaluate(relative, site1, site2, model_params);
            }
          }
        }
      }
    }
  }
  set_energy(energy);
}

}  // namespace feasst
