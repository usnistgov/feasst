#include <vector>
#include "core/include/visit_model.h"
#include "core/include/model_two_body.h"
#include "core/include/model_one_body.h"
#include "core/include/select_list.h"
#include "core/include/constants.h"

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
  DEBUG("visiting model");
  energy_ = 0;
  Position relative;
  double r2 = 0.;
  const Domain& domain = config.domain();
  relative.set_vector(domain.side_length().coord());
  const ModelParams& model_params = config.unique_types().model_params();
  const Select& select_all = config.group_selects()[group_index];
  for (int select1_index = 0;
       select1_index < selection.num_particles();
       ++select1_index) {
    const int part1_index = selection.particle_index(select1_index);
    const Particle part1 = config.select_particle(part1_index);
    TRACE("part1_index " << part1_index << " s " <<
          selection.particle_indices().size() << " " <<
          selection.site_indices().size());
    for (int select2_index = 0;
         select2_index < select_all.num_particles();
         ++select2_index) {
      const int part2_index = select_all.particle_index(select2_index);
      if (part1_index != part2_index) {
        const Particle& part2 = config.select_particle(part2_index);
        for (int site1_index : selection.site_indices(select1_index)) {
          TRACE("site1_index " << site1_index);
          const Site& site1 = part1.sites()[site1_index];
          for (int site2_index : select_all.site_indices(select2_index)) {
            const Site& site2 = part2.sites()[site2_index];
            TRACE("index: " << part1_index << " " << part2_index << " " <<
                  site1_index << " " << site2_index);
            domain.wrap_opt(site1.position(), site2.position(), &relative, &r2);
//    const double squared_distance = relative.squared_distance();
//    const int type1 = site1.type();
//    const int type2 = site2.type();
//    const double cutoff = (*model_params.mixed_cutoff())[type1][type2];
//    TRACE("squared dist: " << squared_distance);
//    if (squared_distance <= cutoff*cutoff) {
            energy_ += model.evaluate(relative, site1, site2, model_params);
//    }
          }
        }
      }
    }
  }
}

void VisitModel::check_energy(const Configuration& config,
    const Model& model,
    const int group_index) {
  model.compute(*this, config, group_index);
  const double en_group = energy();

  // select each particle and compare half the sum with the whole
  SelectList select;
  double en_select = 0;
  const int num = config.num_particles(group_index);
  for (int part = 0; part < num; ++part) {
    select.particle(part, config, group_index);
    model.compute(*this, config, select, group_index);
    TRACE("part " << part << " en " << energy());
    en_select += 0.5*energy();
  }
  ASSERT(std::abs(en_group - en_select) < num*num*1e-15, "Error with " <<
    "visitor implementation. The energy of " <<
    MAX_PRECISION << "group(" << group_index << "): " << en_group << " "
    "is not consistent with half the sum of the energies of the selected " <<
    "particles: " << en_select << ". The difference is: " <<
    en_group - en_select << " with tolerance: " << num*num*1e-15);
}

}  // namespace feasst
