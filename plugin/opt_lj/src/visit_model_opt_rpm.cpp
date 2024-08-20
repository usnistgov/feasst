#include <cmath>  // rint
#include "utils/include/utils.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/model_params.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "configuration/include/physical_constants.h"
#include "system/include/visit_model_inner.h"
#include "opt_lj/include/visit_model_opt_rpm.h"

namespace feasst {

void VisitModelOptRPM::precompute(Configuration * config) {
  VisitModel::precompute(config);
  const ModelParams& existing = config->model_params();
  alpha_ = existing.property("alpha");
  conversion_factor_ = existing.constants().charge_conversion();
  //init_erfc_(existing.cutoff().mixed_max());
}

void VisitModelOptRPM::compute(
    ModelTwoBody * model,
    const ModelParams& model_params,
    const Select& selection,
    Configuration * config,
    const int group_index) {
  FATAL("Don't use VisitModelOptRPM, it has some issues still.");
  DEBUG("visiting model");
  zero_energy();
  double energy = 0.;
  const Domain& domain = config->domain();
  double xi, yi, zi, dx, dy, dz;
  const double lx = domain.side_length(0),
               ly = domain.side_length(1),
               lz = domain.side_length(2);
  const Select& select_all = config->group_select(group_index);
  bool is_old_config = false;
  if (selection.trial_state() == 0 ||
      selection.trial_state() == 2) {
    is_old_config = true;
  }
  Position pbc(3);
  if (selection.num_particles() == 1) {
    for (int select1_index = 0;
         select1_index < selection.num_particles();
         ++select1_index) {
      const int part1_index = selection.particle_index(select1_index);
      const Particle& part1 = config->select_particle(part1_index);
      TRACE("part1_index " << part1_index << " s " <<
            selection.particle_indices().size() << " " <<
            selection.site_indices().size());
      for (int select2_index = 0;
           select2_index < select_all.num_particles();
           ++select2_index) {
        const int part2_index = select_all.particle_index(select2_index);
        const Particle& part2 = config->select_particle(part2_index);
        if (part1_index != part2_index) {
          for (int site1_index : selection.site_indices(select1_index)) {
            TRACE("site1_index " << site1_index);
            const Site& site1 = part1.site(site1_index);
            const int type1 = site1.type();
            const std::vector<double>& coord1 = site1.position().coord();
            xi = coord1[0];
            yi = coord1[1];
            zi = coord1[2];
            for (int site2_index : select_all.site_indices(select2_index)) {
              TRACE("index: " << part1_index << " " << part2_index << " " <<
                    site1_index << " " << site2_index);
              const Site& site2 = part2.site(site2_index);
              if (!is_old_config) get_inner_()->clear_ixn(part1_index, site1_index, part2_index, site2_index);
              const int type2 = site2.type();
              const std::vector<double>& coord2 = site2.position().coord();
              dx = xi - coord2[0];
              dx -= lx*std::rint(dx/lx);
              pbc.set_coord(0, dx);
              dy = yi - coord2[1];
              dy -= ly*std::rint(dy/ly);
              pbc.set_coord(1, dy);
              dz = zi - coord2[2];
              dz -= lz*std::rint(dz/lz);
              pbc.set_coord(2, dz);
              const double squared_distance = dx*dx + dy*dy + dz*dz;
              const double cutoff = model_params.select(cutoff_index()).mixed_values()[type1][type2];
              if (squared_distance <= cutoff*cutoff) {
                const double sigma = model_params.select(sigma_index()).mixed_values()[type1][type2];
                if (squared_distance <= sigma*sigma) {
                  set_energy(NEAR_INFINITY);
                  return;
                }
                const double mixed_charge = model_params.select(charge_index()).mixed_values()[type1][type2];
                const double distance = std::sqrt(squared_distance);
                if (std::abs(distance) < NEAR_ZERO) {
                  set_energy(NEAR_INFINITY);
                  return;
                }
                const double en = mixed_charge*conversion_factor_*std::erfc(alpha_*distance)/distance;
        get_inner_()->update_ixn(en, part1_index, site1_index, type1, part2_index,
                   site2_index, type2, squared_distance, &pbc, is_old_config, *config);
                energy += en;
              }
            }
          }
        }
      }
    }
  } else {
    for (int select2_index = 0;
         select2_index < select_all.num_particles();
         ++select2_index) {
      const int part2_index = select_all.particle_index(select2_index);
      const Particle& part2 = config->select_particle(part2_index);
      if (!find_in_list(part2_index, selection.particle_indices())) {
        for (int select1_index = 0;
             select1_index < selection.num_particles();
             ++select1_index) {
          const int part1_index = selection.particle_index(select1_index);
          const Particle& part1 = config->select_particle(part1_index);
          for (const int site1_index : selection.site_indices(select1_index)) {
            const Site& site1 = part1.site(site1_index);
            const int type1 = site1.type();
            const std::vector<double>& coord1 = site1.position().coord();
            xi = coord1[0];
            yi = coord1[1];
            zi = coord1[2];
            for (const int site2_index : select_all.site_indices(select2_index)) {
              const Site& site2 = part2.site(site2_index);
              if (!is_old_config) get_inner_()->clear_ixn(part1_index, site1_index, part2_index, site2_index);
              const int type2 = site2.type();
              const std::vector<double>& coord2 = site2.position().coord();
              dx = xi - coord2[0];
              dx -= lx*std::rint(dx/lx);
              pbc.set_coord(0, dx);
              dy = yi - coord2[1];
              dy -= ly*std::rint(dy/ly);
              pbc.set_coord(1, dy);
              dz = zi - coord2[2];
              dz -= lz*std::rint(dz/lz);
              pbc.set_coord(2, dz);
              const double squared_distance = dx*dx + dy*dy + dz*dz;
              const double cutoff = model_params.select(cutoff_index()).mixed_values()[type1][type2];
              if (squared_distance <= cutoff*cutoff) {
                const double sigma = model_params.select(sigma_index()).mixed_values()[type1][type2];
                if (squared_distance <= sigma*sigma) {
                  set_energy(NEAR_INFINITY);
                  return;
                }
                const double mixed_charge = model_params.select(charge_index()).mixed_values()[type1][type2];
                const double distance = std::sqrt(squared_distance);
                if (std::abs(distance) < NEAR_ZERO) {
                  set_energy(NEAR_INFINITY);
                  return;
                }
                const double en = mixed_charge*conversion_factor_*std::erfc(alpha_*distance)/distance;
        get_inner_()->update_ixn(en, part1_index, site1_index, type1, part2_index,
                   site2_index, type2, squared_distance, &pbc, is_old_config, *config);
                energy += en;
              }
            }
          }
        }
      }
    }
    for (int select1_index = 0;
         select1_index < selection.num_particles() - 1;
         ++select1_index) {
      const int part1_index = selection.particle_index(select1_index);
      const Particle& part1 = config->select_particle(part1_index);
      for (int select2_index = select1_index + 1;
           select2_index < selection.num_particles();
           ++select2_index) {
        const int part2_index = selection.particle_index(select2_index);
        const Particle& part2 = config->select_particle(part2_index);
        if (part1_index != part2_index) {
          for (const int site1_index : selection.site_indices(select1_index)) {
            const Site& site1 = part1.site(site1_index);
            const int type1 = site1.type();
            const std::vector<double>& coord1 = site1.position().coord();
            xi = coord1[0];
            yi = coord1[1];
            zi = coord1[2];
            for (const int site2_index : selection.site_indices(select2_index)) {
              const Site& site2 = part2.site(site2_index);
              if (!is_old_config) get_inner_()->clear_ixn(part1_index, site1_index, part2_index, site2_index);
              const int type2 = site2.type();
              const std::vector<double>& coord2 = site2.position().coord();
              dx = xi - coord2[0];
              dx -= lx*std::rint(dx/lx);
              pbc.set_coord(0, dx);
              dy = yi - coord2[1];
              dy -= ly*std::rint(dy/ly);
              pbc.set_coord(1, dy);
              dz = zi - coord2[2];
              dz -= lz*std::rint(dz/lz);
              pbc.set_coord(2, dz);
              const double squared_distance = dx*dx + dy*dy + dz*dz;
              const double cutoff = model_params.select(cutoff_index()).mixed_values()[type1][type2];
              if (squared_distance <= cutoff*cutoff) {
                const double sigma = model_params.select(sigma_index()).mixed_values()[type1][type2];
                if (squared_distance <= sigma*sigma) {
                  set_energy(NEAR_INFINITY);
                  return;
                }
                const double mixed_charge = model_params.select(charge_index()).mixed_values()[type1][type2];
                const double distance = std::sqrt(squared_distance);
                if (std::abs(distance) < NEAR_ZERO) {
                  set_energy(NEAR_INFINITY);
                  return;
                }
                const double en = mixed_charge*conversion_factor_*std::erfc(alpha_*distance)/distance;
        get_inner_()->update_ixn(en, part1_index, site1_index, type1, part2_index,
                   site2_index, type2, squared_distance, &pbc, is_old_config, *config);
                energy += en;
              }
            }
          }
        }
      }
    }
  }
  set_energy(energy);
}

class MapVisitModelOptRPM {
 public:
  MapVisitModelOptRPM() {
    VisitModelOptRPM().deserialize_map()["VisitModelOptRPM"] =
      std::make_shared<VisitModelOptRPM>();
  }
};

static MapVisitModelOptRPM mapper_ = MapVisitModelOptRPM();

std::shared_ptr<VisitModel> VisitModelOptRPM::create(std::istream& istr) const {
  return std::make_shared<VisitModelOptRPM>(istr);
}

VisitModelOptRPM::VisitModelOptRPM(std::istream& istr) : VisitModel(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(2047 == version, version);
}

void VisitModelOptRPM::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_visit_model_(ostr);
  feasst_serialize_version(2047, ostr);
}

}  // namespace feasst
