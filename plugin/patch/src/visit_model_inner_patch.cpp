#include "system/include/model_two_body.h"
#include "patch/include/visit_model_inner_patch.h"

namespace feasst {

void VisitModelInnerPatch::compute(
    const int part1_index,
    const int site1_index,
    const int part2_index,
    const int site2_index,
    const Configuration * config,
    const ModelParams& model_params,
    const ModelTwoBody& model,
    const bool is_old_config,
    Position * relative,
    Position * pbc) {
  // HWH copy-pasted from base class.. make this a function?
  if (is_old_config && is_energy_map_queryable()) {
    DEBUG("using old map");
    query_ixn(part1_index, site1_index, part2_index, site2_index);
    return;
  }
  const Particle& part1 = config->select_particle(part1_index);
  const Site& site1 = part1.site(site1_index);
  const Particle& part2 = config->select_particle(part2_index);
  const Site& site2 = part2.site(site2_index);
  clear_ixn(part1_index, site1_index, part2_index, site2_index);
  double squared_distance;
  config->domain()->wrap_opt(site1.position(), site2.position(), relative, pbc, &squared_distance);
  const int type1 = site1.type();
  const int type2 = site2.type();
  const double cutoff = model_params.mixed_cutoff()[type1][type2];
  // this cutoff should be based on max possible between centers.
  // HWH add error check for this assumption
//  TRACE("mode parm ang " << model_params.select("angle")->mixed_value(1, 1));
  if (squared_distance <= cutoff*cutoff) {
    // loop through bonds on center1
    const int part1_type = part1.type();
    for (const int dir1_index : config->particle_type(part1_type).bond_neighbor()[site1_index]) {
      const Site& dir1 = part1.site(dir1_index);
      if (dir1.is_physical()) {
        const int dir1_type = dir1.type();
        if (config->unique_type(part1_type).site(dir1_type).is_director()) {
          // loop through bonds on center2
          const int part2_type = part2.type();
          for (const int dir2_index : config->particle_type(part2_type).bond_neighbor()[site2_index]) {
            const Site& dir2 = part2.site(dir2_index);
            if (dir2.is_physical()) {
              const int dir2_type = dir2.type();
              if (config->unique_type(part2_type).site(dir2_type).is_director()) {
                const double dircut = model_params.mixed_cutoff()[dir1_type][dir2_type];
                if (squared_distance <= dircut*dircut) {
                  Position dir1_pos = dir1.position();
                  dir1_pos.subtract(site1.position());
                  dir1_pos.multiply(-1.);
                  const double dir1_sq_length = dir1_pos.squared_distance();
                  const double cosp1 = dir1_pos.dot_product(*relative)/sqrt(squared_distance*dir1_sq_length);
                  TRACE("cosp1 " << cosp1 << " cosacut " << cos_patch_angle_.value(dir1_type));
                  if (cosp1 >= cos_patch_angle_.value(dir1_type)) {
                    Position dir2_pos = dir2.position();
                    dir2_pos.subtract(site2.position());
                    const double dir2_sq_length = dir2_pos.squared_distance();
                    const double cosp2 = dir2_pos.dot_product(*relative)/sqrt(squared_distance*dir2_sq_length);
                    TRACE("cosp2 " << cosp2 << " cosacut " << cos_patch_angle_.value(dir2_type));
                    if (cosp2 >= cos_patch_angle_.value(dir2_type)) {
                      const double en = model.energy(squared_distance, dir1_type, dir2_type, model_params);
                      update_ixn(en, part1_index, site1_index, part2_index,
                                 site2_index, squared_distance, pbc);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

class MapVisitModelInnerPatch {
 public:
  MapVisitModelInnerPatch() {
    VisitModelInnerPatch().deserialize_map()["VisitModelInnerPatch"] =
      std::make_shared<VisitModelInnerPatch>();
  }
};

static MapVisitModelInnerPatch mapper_ = MapVisitModelInnerPatch();

VisitModelInnerPatch::VisitModelInnerPatch(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 255, "unrecognized version: " << version);
  feasst_deserialize_fstobj(&cos_patch_angle_, istr);
}

void VisitModelInnerPatch::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  feasst_serialize_version(255, ostr);
  feasst_serialize_fstobj(cos_patch_angle_, ostr);
}

}  // namespace feasst
