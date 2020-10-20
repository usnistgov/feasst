#include <cmath>  // rint
#include "utils/include/serialize.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "system/include/lennard_jones.h"
#include "opt_lj/include/visit_model_opt_lj.h"

namespace feasst {

void VisitModelOptLJ::compute(
    const ModelTwoBody& model,
    const ModelParams& model_params,
    const Select& selection,
    Configuration * config,
    const int group_index) {
  DEBUG("visiting model");
  zero_energy();
  double energy = 0.;
  const Domain& domain = config->domain();
  double xi, yi, zi, dx, dy, dz;
  const double lx = domain.side_length(0),
               ly = domain.side_length(1),
               lz = domain.side_length(2);
  const Select& select_all = config->group_selects()[group_index];
//  LennardJones lj_; // model appears to be faster on stack
  // HWH implement multi-particle selection by sorting group selection
  // for particles that are in both selectiona nd group_index.
  // treat those particles separately so no double counting.
  // then remove the part1 != part2 check
  ASSERT(selection.num_particles() == 1, "for multiparticle selections " <<
    "implement a separate loop for particles in both group and selection. " <<
    "Select: " << selection.str());
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
            const int type2 = site2.type();
            const std::vector<double>& coord2 = site2.position().coord();
            dx = xi - coord2[0];
            dx -= lx*std::rint(dx/lx);
            dy = yi - coord2[1];
            dy -= ly*std::rint(dy/ly);
            dz = zi - coord2[2];
            dz -= lz*std::rint(dz/lz);
            const double squared_distance = dx*dx + dy*dy + dz*dz;
            const double cutoff = model_params.mixed_cutoff()[type1][type2];
            if (squared_distance <= cutoff*cutoff) {
//                  const double en = lj_.energy(squared_distance,
//                    type1,
//                    type2,
//                    model_params);
              const double sigma = model_params.mixed_sigma()[type1][type2];
              const double sigma_squared = sigma*sigma;
              const double epsilon = model_params.mixed_epsilon()[type1][type2];
              const double rinv2 = sigma_squared/squared_distance;
              const double rinv6 = rinv2*rinv2*rinv2;
              const double en = 4.*epsilon*rinv6*(rinv6 - 1.);
              energy += en;
            }
          }
        }
      }
    }
  }
  set_energy(energy);
}

class MapVisitModelOptLJ {
 public:
  MapVisitModelOptLJ() {
    VisitModelOptLJ().deserialize_map()["VisitModelOptLJ"] =
      std::make_shared<VisitModelOptLJ>();
  }
};

static MapVisitModelOptLJ mapper_ = MapVisitModelOptLJ();

std::shared_ptr<VisitModel> VisitModelOptLJ::create(std::istream& istr) const {
  return std::make_shared<VisitModelOptLJ>(istr);
}

VisitModelOptLJ::VisitModelOptLJ(std::istream& istr) : VisitModel(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(825 == version, version);
}

void VisitModelOptLJ::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_visit_model_(ostr);
  feasst_serialize_version(825, ostr);
}

}  // namespace feasst
