#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/model_params.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "system/include/model_two_body.h"
#include "patch/include/visit_model_inner_patch.h"

namespace feasst {

VisitModelInnerPatch::VisitModelInnerPatch(argtype * args) : VisitModelInner(args) {
  class_name_ = "VisitModelInnerPatch";
}
VisitModelInnerPatch::VisitModelInnerPatch(argtype args) : VisitModelInnerPatch(&args) {
  feasst_check_all_used(args);
}

void VisitModelInnerPatch::compute(
    const int part1_index,
    const int site1_index,
    const int part2_index,
    const int site2_index,
    const Configuration * config,
    const ModelParams& model_params,
    ModelTwoBody * model,
    const bool is_old_config,
    Position * relative,
    Position * pbc,
    const double weight) {
  TRACE("part1_index " << part1_index);
  TRACE("part2_index " << part2_index);
  TRACE("site1_index " << site1_index);
  TRACE("site2_index " << site2_index);
  const Particle& part1 = config->select_particle(part1_index);
  const Site& site1 = part1.site(site1_index);
  const Particle& part2 = config->select_particle(part2_index);
  const Site& site2 = part2.site(site2_index);
  clear_ixn(part1_index, site1_index, part2_index, site2_index);
  const ModelParam& director = model_params.select(director_index_);
//  const ModelParam& cos_patch_angle = model_params.select(cos_patch_angle_index_);
  double squared_distance;
  config->domain().wrap_opt(site1.position(), site2.position(), relative, pbc, &squared_distance);
  const int type1 = site1.type();
  const int type2 = site2.type();
  const double cutoff = model_params.select(cutoff_index()).mixed_values()[type1][type2];
  // this cutoff should be based on max possible between centers.
  // HWH add error check for this assumption
  //TRACE("mode parm ang " << model_params.select("patch_angle").mixed_value(1, 1));
  TRACE("cutoff " << cutoff);
  TRACE("squared_distance " << squared_distance);
  if (squared_distance <= cutoff*cutoff &&
      squared_distance > 0.) {
    // loop through bonds on center1
    const int part1_type = part1.type();
    for (const int dir1_index : config->particle_type(part1_type).bond_neighbors(site1_index)) {
      TRACE("dir1_index " << dir1_index);
      const Site& dir1 = part1.site(dir1_index);
      if (dir1.is_physical()) {
        const int dir1_type = dir1.type();
        if (director.value(dir1_type) > 0.5) {
          // loop through bonds on center2
          const int part2_type = part2.type();
          for (const int dir2_index : config->particle_type(part2_type).bond_neighbors(site2_index)) {
            TRACE("dir2_index " << dir2_index);
            const Site& dir2 = part2.site(dir2_index);
            if (dir2.is_physical()) {
              const int dir2_type = dir2.type();
              TRACE("dir2_type " << dir2_type);
              if (director.value(dir2_type) > 0.5) {
                const double dircut = model_params.select(cutoff_index()).mixed_values()[dir1_type][dir2_type];
                TRACE("dircut " << dircut);
                if (squared_distance <= dircut*dircut) {
                  dir1_pos_.set_vector(dir1.position().coord());
                  dir1_pos_.subtract(site1.position());
                  dir1_pos_.multiply(-1.);
                  const double dir1_sq_length = dir1_pos_.squared_distance();
                  TRACE("dir1_sq_length " << dir1_sq_length);
                  TRACE("sqdist " << squared_distance);
                  const double cosp1 = dir1_pos_.dot_product(*relative)/std::sqrt(squared_distance*dir1_sq_length);
                  TRACE("cosp1 " << cosp1 << " cosacut " << cos_patch_angle_.value(dir1_type));
                  if (cosp1 >= cos_patch_angle_.value(dir1_type)) {
                    dir2_pos_.set_vector(dir2.position().coord());
                    dir2_pos_.subtract(site2.position());
                    const double dir2_sq_length = dir2_pos_.squared_distance();
                    const double cosp2 = dir2_pos_.dot_product(*relative)/std::sqrt(squared_distance*dir2_sq_length);
                    TRACE("cosp2 " << cosp2 << " cosacut " << cos_patch_angle_.value(dir2_type));
                    if (cosp2 >= cos_patch_angle_.value(dir2_type)) {
                      const double en = weight*model->energy(squared_distance, dir1_type, dir2_type, model_params);
                      TRACE("en " << en);
                      update_ixn(en, part1_index, site1_index, type1, part2_index,
                                 site2_index, type2, squared_distance, pbc, is_old_config, *config);
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

FEASST_MAPPER(VisitModelInnerPatch,);

VisitModelInnerPatch::VisitModelInnerPatch(std::istream& istr) : VisitModelInner(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 255, "unrecognized version: " << version);
  feasst_deserialize_fstobj(&cos_patch_angle_, istr);
  feasst_deserialize(&director_index_, istr);
  feasst_deserialize_fstobj(&dir1_pos_, istr);
  feasst_deserialize_fstobj(&dir2_pos_, istr);
}

void VisitModelInnerPatch::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_visit_model_inner_(ostr);
  feasst_serialize_version(255, ostr);
  feasst_serialize_fstobj(cos_patch_angle_, ostr);
  feasst_serialize(director_index_, ostr);
  feasst_serialize_fstobj(dir1_pos_, ostr);
  feasst_serialize_fstobj(dir2_pos_, ostr);
}

void VisitModelInnerPatch::precompute(Configuration * config) {
  VisitModelInner::precompute(config);
  director_index_ = config->model_params().index("director");
  //config->add(std::make_shared<CosPatchAngle>());
  cos_patch_angle_.set_param(config->model_params());
  for (int type = 0; type < config->num_site_types(); ++type) {
    cos_patch_angle_.compute(type, config->model_params());
  }
  //cos_patch_angle_index_ = config->model_params().index("cos_patch_angle");
  //INFO("cos_patch_angle_index_ " << cos_patch_angle_index_);
  dir1_pos_.set_to_origin(config->domain().dimension());
  dir2_pos_.set_to_origin(config->domain().dimension());
}

}  // namespace feasst
