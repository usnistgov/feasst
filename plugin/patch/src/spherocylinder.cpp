#include <cmath>
#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "system/include/model_two_body.h"
#include "patch/include/spherocylinder.h"

namespace feasst {

Spherocylinder::Spherocylinder(argtype * args) : VisitModelInner(args) {
  class_name_ = "Spherocylinder";
}
Spherocylinder::Spherocylinder(argtype args) : Spherocylinder(&args) {
  feasst_check_all_used(args);
}

double Spherocylinder::calc_sph_sq_dist_vega_(const double rr, const double rw1,
    const double rw2, const double w1w2, const double lh1, const double lh2) {
  const double cc = 1 - w1w2*w1w2;
  if (cc < 1e-6) {
    TRACE("parallel");
    if (rw1 != 0 && rw2 != 0) {
    //if (rw1 && rw2) {
       xla_ = rw1/2;
       xmu_ = -rw2/2;
      } else {
        return rr;
      }
  } else {
    TRACE("not parallel");
    xla_ = (rw1 - w1w2*rw2)/cc;
    xmu_ = (-rw2 + w1w2*rw1)/cc;
  }

  if (std::abs(xla_) > lh1 || std::abs(xmu_) > lh2) {
    if (std::abs(xla_) - lh1 > std::abs(xmu_) - lh2) {
      xla_= lh1*sgn(xla_);
      xmu_= xla_*w1w2 - rw2;
      if (std::abs(xmu_) > lh2) {
        xmu_= lh2*sgn(xmu_);
      }
    } else {
      xmu_= lh2*sgn(xmu_);
      xla_= xmu_*w1w2 + rw1;
      if (std::abs(xla_) > lh1) {
        xla_ = lh1*sgn(xla_);
      }
    }
  }
  return rr + xla_*xla_ + xmu_*xmu_ + 2*(xmu_*rw2 - xla_*(rw1 + xmu_*w1w2));
}

void Spherocylinder::compute(
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
  set_interacted(0);
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
  double squared_distance;
  TRACE("site1pos " << site1.position().str());
  TRACE("site2pos " << site2.position().str());
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
      squared_distance >= 0.) {
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
                dir1_pos_.set_vector(dir1.position().coord());
                dir1_pos_.subtract(site1.position());
                dir1_pos_.normalize();
                TRACE("r " << relative->str());
                TRACE("w1 " << dir1_pos_.str());
                TRACE("sqdist " << squared_distance);
                dir2_pos_.set_vector(dir2.position().coord());
                dir2_pos_.subtract(site2.position());
                dir2_pos_.normalize();
                const double lh1 = 0.5*length_.value(dir1_type);
                TRACE("w2 " << dir2_pos_.str());
                const double lh2 = 0.5*length_.value(dir2_type);
                const double rr = squared_distance;
                TRACE("rr " << rr);
                const double rw1 = relative->dot_product(dir1_pos_);
                TRACE("rw1 " << rw1);
                const double rw2 = relative->dot_product(dir2_pos_);
                TRACE("rw2 " << rw2);
                const double w1w2 = dir1_pos_.dot_product(dir2_pos_);
                TRACE("w1w2 " << w1w2);
                const double sph_sq_dist = calc_sph_sq_dist_vega_(rr, rw1, rw2, w1w2, lh1, lh2);
                TRACE("sph_sq_dist " << sph_sq_dist);
                TRACE("sph_dist " << std::sqrt(sph_sq_dist));
                const double dircut = model_params.select(cutoff_index()).mixed_values()[dir1_type][dir2_type];
                TRACE("dircut " << dircut);
                if (sph_sq_dist <= dircut*dircut) {
                  const double en = weight*model->energy(sph_sq_dist, dir1_type, dir2_type, model_params);
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

FEASST_MAPPER(Spherocylinder,);

Spherocylinder::Spherocylinder(std::istream& istr) : VisitModelInner(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 6768, "unrecognized version: " << version);
  feasst_deserialize_fstobj(&length_, istr);
  feasst_deserialize(&director_index_, istr);
  feasst_deserialize_fstobj(&dir1_pos_, istr);
  feasst_deserialize_fstobj(&dir2_pos_, istr);
}

void Spherocylinder::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_visit_model_inner_(ostr);
  feasst_serialize_version(6768, ostr);
  feasst_serialize_fstobj(length_, ostr);
  feasst_serialize(director_index_, ostr);
  feasst_serialize_fstobj(dir1_pos_, ostr);
  feasst_serialize_fstobj(dir2_pos_, ostr);
}

void Spherocylinder::precompute(Configuration * config) {
  VisitModelInner::precompute(config);
  director_index_ = config->model_params().index("director");
  length_.set_param(config->model_params());
  for (int type = 0; type < config->num_site_types(); ++type) {
    length_.compute(type, config->model_params());
  }
  dir1_pos_.set_to_origin(config->domain().dimension());
  dir2_pos_.set_to_origin(config->domain().dimension());
}

}  // namespace feasst
