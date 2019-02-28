#include "patch/include/visit_model_patch.h"
#include "core/include/model_two_body.h"
#include "core/include/constants.h"

namespace feasst {

void VisitModelPatch::inner_(
    const int part1_index,
    const int site1_index,
    const int part2_index,
    const int site2_index,
    const Configuration * config,
    const ModelParams& model_params,
    const ModelTwoBody& model,
    Position * relative) {
  const Particle& part1 = config->select_particle(part1_index);
  const Site& site1 = part1.site(site1_index);
  const Particle& part2 = config->select_particle(part2_index);
  const Site& site2 = part2.site(site2_index);
  double squared_distance;
  config->domain().wrap_opt(site1.position(), site2.position(), relative, &squared_distance);
  const int type1 = site1.type();
  const int type2 = site2.type();
  const double cutoff = model_params.mixed_cutoff()[type1][type2];
  // this cutoff should be based on max possible between centers.
  // HWH add error check for this assumption
  if (squared_distance <= cutoff*cutoff) {
    // loop through bonds on center1
    const int part1_type = part1.type();
    for (const int dir1_index : config->particle_type(part1_type).bond_neighbor()[site1_index]) {
      const Site& dir1 = part1.site(dir1_index);
      const int dir1_type = dir1.type();
      if (config->unique_type(part1_type).site(dir1_type).is_director()) {
        // loop through bonds on center2
        const int part2_type = part2.type();
        for (const int dir2_index : config->particle_type(part2_type).bond_neighbor()[site2_index]) {
          const Site& dir2 = part2.site(dir2_index);
          const int dir2_type = dir2.type();
          if (config->unique_type(part2_type).site(dir2_type).is_director()) {
            const double dircut = model_params.mixed_cutoff()[dir1_type][dir2_type];
            if (squared_distance <= dircut*dircut) {
              Position dir1_pos = dir1.position();
              dir1_pos.subtract(site1.position());
              dir1_pos.multiply(-1.);
              const double dir1_sq_length = dir1_pos.squared_distance();
              const double cosp1 = dir1_pos.dot_product(*relative);
              TRACE("cosp1 " << cosp1 << " cpa_sq_ " << cpa_sq_ << " cosa2 " << cosp1*cosp1/squared_distance/dir1_sq_length);
              if (cosp1 >= 0 && cosp1*cosp1/squared_distance/dir1_sq_length >= cpa_sq_) {
                Position dir2_pos = dir2.position();
                dir2_pos.subtract(site2.position());
                const double dir2_sq_length = dir2_pos.squared_distance();
                const double cosp2 = dir2_pos.dot_product(*relative);
                TRACE("cosp2 " << cosp2 << " cpa_sq_ " << cpa_sq_ << " cosa2 " << cosp2*cosp2/squared_distance/dir2_sq_length);
                if (cosp2 >= 0 && cosp2*cosp2/squared_distance/dir2_sq_length >= cpa_sq_) {
                  const double en = model.energy(squared_distance, dir1_type, dir2_type, model_params);
                  increment_energy(en);
                }
              }
            }
          }
        }
      }
    }
  }
}

}  // namespace feasst
