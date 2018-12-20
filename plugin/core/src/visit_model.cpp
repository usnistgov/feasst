#include <vector>
#include "core/include/visit_model_cell.h"
#include "core/include/model_two_body.h"
#include "core/include/model_one_body.h"

namespace feasst {

void VisitModel::compute(const Configuration& config,
    const ModelOneBody& model,
    const int group_index) {
  energy_ = 0.;
  const ModelParams& model_params = config.unique_types().model_params();
  for (int part_index = 0;
       part_index < config.num_particles(group_index);
       ++part_index) {
    const Particle& ppart = config.particle(part_index);
    for (const Site& site : ppart.sites()) {
      energy_ += model.evaluate(site, config, model_params);
    }
  }
}

void VisitModel::compute(const Configuration& config,
                         const ModelOneBody& model,
                         const Select& selection,
                         const int group_index) {
  ASSERT(group_index == 0, "not implemented because redundant to selection");
  energy_ = 0;
  ERROR("HWH: add wrapping of site positions");
  const ModelParams& model_params = config.unique_types().model_params();
  for (int select_index = 0; select_index < selection.num_particles(); ++select_index) {
    const int particle_index = selection.particle_index(select_index);
    const Particle& part = config.select_particle(particle_index);
    for (int site_index : selection.site_indices(select_index)) {
      const Site& site = part.sites()[site_index];
      energy_ += model.evaluate(site, config, model_params);
    }
  }
}

void VisitModel::compute(const Configuration& config,
    const ModelTwoBody& model,
    const int group_index) {
  energy_ = 0;
  Position relative;
  double r2 = 0.;
  const Domain& domain = config.domain();
  relative.set_vector(domain.side_length().coord());
  const ModelParams& model_params = config.unique_types().model_params();
  for (int part1_index = 0;
       part1_index < config.num_particles(group_index) - 1;
       ++part1_index) {
    const Particle part1 = config.particle(part1_index);
    for (const Site& site1 : part1.sites()) {
      for (int part2_index = part1_index + 1;
           part2_index < config.num_particles(group_index);
           ++part2_index) {
        const Particle part2 = config.particle(part2_index);
        for (const Site& site2 : part2.sites()) {
          domain.wrap_opt(site1.position(), site2.position(), &relative, &r2);
          energy_ += model.evaluate(relative, site1, site2, model_params);
        }
      }
    }
  }
}

void VisitModel::compute(const Configuration& config,
                         const ModelTwoBody& model,
                         const Select& selection,
                         const int group_index) {
  DEBUG("compute");
  energy_ = 0;
  Position relative;
  double r2 = 0.;
  const Domain& domain = config.domain();
  relative.set_vector(domain.side_length().coord());
  const ModelParams& model_params = config.unique_types().model_params();
  const Select& select_all = config.group_selects()[group_index];
  for (int select_index = 0; select_index < selection.num_particles(); ++select_index) {
    const int part1_index = selection.particle_index(select_index);
    const Particle part1 = config.select_particle(part1_index);
    // selection.check_size();
    TRACE("part1_index " << part1_index << " s " << selection.particle_indices().size() << " " << selection.site_indices().size());
    for (int site1_index : selection.site_indices(select_index)) {
      TRACE("site1_index " << site1_index);
      const Site& site1 = part1.sites()[site1_index];
      for (int select2_index = 0; select2_index < select_all.num_particles(); ++select2_index) {
        const int part2_index = select_all.particle_index(select2_index);
        if (part1_index != part2_index) {
          const Particle& part2 = config.select_particle(part2_index);
          for (const Site& site2 : part2.sites()) {
            domain.wrap_opt(site1.position(), site2.position(), &relative, &r2);
            energy_ += model.evaluate(relative, site1, site2, model_params);
          }
        }
      }
    }
  }
}

}  // namespace feasst
