#include <vector>
#include "core/include/visit_model.h"
#include "core/include/model_two_body.h"
#include "core/include/model_one_body.h"

namespace feasst {

void VisitModel::loop_by_particle(const Configuration& config,
                                  const ModelOneBody& model,
                                  const int iPart) {
  energy_ = kloop_by_particle(config, model, iPart);
}

double VisitModel::kloop_by_particle(const Configuration& config,
                                     const ModelOneBody& model,
                                     const int iPart) const {
  double energy = 0.;
  int iPartBegin = iPart, iPartEnd = iPart + 1;
  if (iPart == -1) {
    iPartBegin = 0;
    iPartEnd = config.num_particles();
  }
  const ModelParams& model_params = config.unique_types().model_params();
  for (int part = iPartBegin; part < iPartEnd; ++part) {
    const Particle& ppart = config.particle(part);
    for (const Site& site : ppart.sites()) {
      energy += model.evaluate(site, config, model_params);
    }
  }
  return energy;
}

void VisitModel::energy_of_selection(const Configuration& config,
                                     const ModelOneBody& model) {
  energy_ = 0;
  const ModelParams& model_params = config.unique_types().model_params();
  const Selection& selection = config.selection();
  for (int select_index = 0; select_index < selection.num(); ++select_index) {
    const int particle_index = selection.particle_index(select_index);
    const Particle& part = config.particle(particle_index);
    for (int site_index : config.selection().site_indices(select_index)) {
      const Site& site = part.sites()[site_index];
      energy_ += model.evaluate(site, config, model_params);
    }
  }
}

//  void energy_of_selection(const Configuration& config,
//                           const ModelTwoBody& model);

void VisitModel::loop_by_particle(const Configuration& config,
                                  const ModelTwoBody& model,
                                  const int iPart) {
  energy_ = kloop_by_particle(config, model, iPart);
}

double VisitModel::kloop_by_particle(const Configuration& config,
                                     const ModelTwoBody& model,
                                     const int iPart) const {
  DEBUG("kloop_by_particle");
  double energy = 0.;
  Position difference;
  const DomainCuboid &domain = config.domain();
  const ModelParams& model_params = config.unique_types().model_params();
  int iPartBegin = iPart, iPartEnd = iPart + 1;
  if (iPart == -1) {
    iPartBegin = 0;
    iPartEnd = config.num_particles() - 1;
  }
  for (int part1 = iPartBegin; part1 < iPartEnd; ++part1) {
    const Particle& ppart1 = config.particle(part1);
    for (const Site& site1 : ppart1.sites()) {
      int jPartBegin = 0;
      if (iPart == -1) {
        jPartBegin = part1 + 1;
      }
      for (int part2 = jPartBegin; part2 < config.num_particles(); ++part2) {
        if (part1 != part2) {
          const Particle& ppart2 = config.particle(part2);
          for (const Site& site2 : ppart2.sites()) {
            difference = site1.position();
            difference.subtract(site2.position());
            domain.wrap(&difference);
            TRACE("particle indices: " << part1 << " " << part2);
            energy += model.evaluate(difference, site1, site2, model_params);
          }
        }
      }
    }
  }
  return energy;
}

void VisitModel::energy_of_selection(const Configuration& config,
                                     const ModelTwoBody& model) {
  if (optimization_ == 2) {
    benchmark_(config, model);
    return;
  }
  DEBUG("energy_of_selection");
  energy_ = 0;
  Position relative;
  double r2 = 0.;
  const DomainCuboid& domain = config.domain();
  relative.set_vector(domain.side_length().coord());
  const ModelParams& model_params = config.unique_types().model_params();
  const Selection& selection = config.selection();
  for (int select_index = 0; select_index < selection.num(); ++select_index) {
    const int part1_index = selection.particle_index(select_index);
    const Particle part1 = config.particle(part1_index);
    selection.check_size();
    TRACE("part1_index " << part1_index << " s " << selection.particle_indices().size() << " " << selection.site_indices().size());
    for (int site1_index : config.selection().site_indices(select_index)) {
      TRACE("site1_index " << site1_index);
      const Site& site1 = part1.sites()[site1_index];
      for (int part2_index = 0;
           part2_index < config.num_particles();
           ++part2_index) {
        if (part1_index != part2_index) {
          const Particle& part2 = config.particle(part2_index);
          for (const Site& site2 : part2.sites()) {
            domain.wrap_opt(site1.position(), site2.position(), &relative, &r2);
            energy_ += model.evaluate(relative, site1, site2, model_params);
          }
        }
      }
    }
  }
}

void VisitModel::benchmark_(const Configuration& config,
                            const ModelTwoBody& model) {
  DEBUG("benchmark_");
  ERROR("HWH depreciated");
  energy_ = 0;
  const ModelParams& model_params = config.unique_types().model_params();
  const DomainCuboid &domain = config.domain();
  Position relative;
  double r2;
  const int particle1_index = config.selection().particle_index(0);
  const Site& site1 = config.particle(particle1_index).site(0);
  const std::vector<double>& x1 = site1.position().coord();
  const int dimension = x1.size();
  std::vector<double> dxv(dimension);
  std::vector<double> side_length = domain.side_length().coord();
  relative.set_vector(domain.side_length().coord());
  for (int particle2_index = 0;
       particle2_index < config.num_particles();
       ++particle2_index) {
    if (particle1_index != particle2_index) {
      const Site& site2 = config.particle(particle2_index).site(0);
      domain.wrap_opt(site1.position(), site2.position(), &relative, &r2);
      energy_ += model.evaluate(relative, site1, site2, model_params);
    }
  }
}

}  // namespace feasst
