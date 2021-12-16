#include <cmath>
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "system/include/model_two_body.h"
#include "patch/include/visit_model_inner_patch.h"

namespace feasst {

VisitModelInnerPatch::VisitModelInnerPatch(argtype args) {
  class_name_ = "VisitModelInnerPatch";
  // loop through all arguments, looking for prefix
  const std::string prefix = "patch_degrees_of_type";
  DEBUG("prefix " << prefix);
  DEBUG("prefix sz " << prefix.size());
  argtype copy = args;
  for (std::map<std::string, std::string>::iterator iter = copy.begin();
       iter != copy.end(); ++iter) {
  //for (std::map<std::string, std::string>::reverse_iterator iter = copy.rbegin();
  //     iter != copy.rend(); ++iter) {
    DEBUG("key " << iter->first);
    if (is_found_in(iter->first, prefix)) {
      DEBUG("prefix " << prefix << " found");
      std::string type(iter->first);
      type.erase(type.begin(), type.begin() + prefix.size());
      DEBUG("type " << type);
      const int t = str_to_int(type);
      DEBUG("t " << t);
      const double degrees = dble(iter->first, &args);
      DEBUG("degrees " << degrees);
      const double cosa = std::cos(degrees_to_radians(degrees));
      DEBUG("cosa " << cosa);
      cpatch_override_.push_back({t, cosa});
    }
  }
  check_all_used(args);
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
    Position * pbc) {
//  std::cout << " Info 1 part1_index " << part1_index << std::endl;
//  std::cout << " Info 1 site1_index " << site1_index << std::endl;
//  std::cout << " Info 2 part2_index " << part2_index << std::endl;
//  std::cout << " Info 2 site2_index " << site2_index << std::endl;
  TRACE("part1_index " << part1_index);
  TRACE("part2_index " << part2_index);
  TRACE("site1_index " << site1_index);
  TRACE("site2_index " << site2_index);
  const Particle& part1 = config->select_particle(part1_index);
  const Site& site1 = part1.site(site1_index);
  const Particle& part2 = config->select_particle(part2_index);
  const Site& site2 = part2.site(site2_index);
  clear_ixn(part1_index, site1_index, part2_index, site2_index);
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
        if (director_.value(dir1_type) > 0.5) {
          // loop through bonds on center2
          const int part2_type = part2.type();
          for (const int dir2_index : config->particle_type(part2_type).bond_neighbors(site2_index)) {
            TRACE("dir2_index " << dir2_index);
            const Site& dir2 = part2.site(dir2_index);
            if (dir2.is_physical()) {
              const int dir2_type = dir2.type();
              TRACE("dir2_type " << dir2_type);
              if (director_.value(dir2_type) > 0.5) {
                const double dircut = model_params.select(cutoff_index()).mixed_values()[dir1_type][dir2_type];
                TRACE("dircut " << dircut);
                if (squared_distance <= dircut*dircut) {
                  Position dir1_pos = dir1.position();
                  dir1_pos.subtract(site1.position());
                  dir1_pos.multiply(-1.);
                  const double dir1_sq_length = dir1_pos.squared_distance();
                  TRACE("dir1_sq_length " << dir1_sq_length);
                  TRACE("sqdist " << squared_distance);
                  const double cosp1 = dir1_pos.dot_product(*relative)/std::sqrt(squared_distance*dir1_sq_length);
                  TRACE("cosp1 " << cosp1 << " cosacut " << cos_patch_angle_.value(dir1_type));
                  if (cosp1 >= cos_patch_angle_.value(dir1_type)) {
                    Position dir2_pos = dir2.position();
                    dir2_pos.subtract(site2.position());
                    const double dir2_sq_length = dir2_pos.squared_distance();
                    const double cosp2 = dir2_pos.dot_product(*relative)/std::sqrt(squared_distance*dir2_sq_length);
                    TRACE("cosp2 " << cosp2 << " cosacut " << cos_patch_angle_.value(dir2_type));
                    if (cosp2 >= cos_patch_angle_.value(dir2_type)) {
                      const double en = model->energy(squared_distance, dir1_type, dir2_type, model_params);
                      TRACE("en " << en);
                      update_ixn(en, part1_index, site1_index, type1, part2_index,
                                 site2_index, type2, squared_distance, pbc, is_old_config);
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

VisitModelInnerPatch::VisitModelInnerPatch(std::istream& istr) : VisitModelInner(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 255, "unrecognized version: " << version);
  feasst_deserialize_fstobj(&cos_patch_angle_, istr);
  feasst_deserialize_fstobj(&director_, istr);
  feasst_deserialize(&cpatch_override_, istr);
}

void VisitModelInnerPatch::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_visit_model_inner_(ostr);
  feasst_serialize_version(255, ostr);
  feasst_serialize_fstobj(cos_patch_angle_, ostr);
  feasst_serialize_fstobj(director_, ostr);
  feasst_serialize(cpatch_override_, ostr);
}

void VisitModelInnerPatch::precompute(Configuration * config) {
  VisitModelInner::precompute(config);
  config->add(std::make_shared<PatchAngle>());
  config->add(std::make_shared<CosPatchAngle>());
  config->add(std::make_shared<Director>());
  DEBUG("setting params");
  cos_patch_angle_.set_param(config->model_params());
  director_.set_param(config->model_params());
  for (const std::pair<int, double>& patch : cpatch_override_) {
    DEBUG("patch " << patch.first << " " << patch.second);
    cos_patch_angle_.set(patch.first, patch.second);
  }
  DEBUG(cos_patch_angle_.str());
  DEBUG(director_.str());
}

//void VisitModelInnerPatch::set_patch_angle(const int type,
//    const double degrees) {
//  const double cosa = std::cos(degrees_to_radians(degrees));
//  cos_patch_angle_.set(type, cosa);
//}

}  // namespace feasst
