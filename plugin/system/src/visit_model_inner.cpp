#include <vector>
#include <cmath>
#include "utils/include/serialize.h"
#include "configuration/include/select.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "system/include/visit_model_inner.h"
#include "system/include/model_two_body.h"

namespace feasst {

class MapVisitModelInner {
 public:
  MapVisitModelInner() {
    VisitModelInner().deserialize_map()["VisitModelInner"] =
      std::make_shared<VisitModelInner>();
  }
};

static MapVisitModelInner mapper_visit_model_inner_ = MapVisitModelInner();

std::map<std::string, std::shared_ptr<VisitModelInner> >&
    VisitModelInner::deserialize_map() {
  static std::map<std::string, std::shared_ptr<VisitModelInner> >* ans =
     new std::map<std::string, std::shared_ptr<VisitModelInner> >();
  return *ans;
}

void VisitModelInner::compute(
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
  TRACE("VisitModelInner");
//  if (is_old_config && is_energy_map_queryable()) {
//    DEBUG("using old map");
//    query_ixn(part1_index, site1_index, part2_index, site2_index);
//    return;
//  }
  const Particle& part1 = config->select_particle(part1_index);
  const Site& site1 = part1.site(site1_index);
  if (!is_old_config) clear_ixn(part1_index, site1_index, part2_index, site2_index);
  if (site1.is_physical()) {
    const Particle& part2 = config->select_particle(part2_index);
    const Site& site2 = part2.site(site2_index);
    if (site2.is_physical()) {
      config->domain().wrap_opt(site1.position(), site2.position(), relative,
                                pbc, &squared_distance_);
      const int type1 = site1.type();
      const int type2 = site2.type();
      const double cutoff = model_params.select(cutoff_index()).mixed_values()[type1][type2];
      TRACE("cutoff " << cutoff);
      TRACE("squared_distance_ " << squared_distance_);
      TRACE("indices " << part1_index << " " << site1_index << " " <<
          part2_index << " " << site2_index);
      if (squared_distance_ <= cutoff*cutoff) {
        const double energy = model->energy(squared_distance_, type1, type2,
                                           model_params);
        update_ixn(energy, part1_index, site1_index, type1, part2_index,
                   site2_index, type2, squared_distance_, pbc, is_old_config);
        TRACE("energy " << energy_ << " " << energy);
      } else {
        // if distance is greater than cutoff+outer, then skip the entire
        // particle.
//        TRACE("cutoff_outer_inter " << cutoff_outer_index_);
        if (cutoff_outer_index_ != -1) {
          if (site1_index == 0 && site2_index == 0) {
            const double outer = model_params.select(cutoff_outer_index_).mixed_values()[type1][type2];
            //TRACE("outer " << outer);
            if (outer > 0) {
              if (squared_distance_ > std::pow(cutoff+2.*outer, 2)) {
                TRACE("skipping! dist: " << squared_distance_ << " > " << std::pow(cutoff+outer, 2));
                skip_particle_ = true;
              }
            }
          }
        }
      }
    }
  }
}

bool VisitModelInner::is_energy_map_queryable() const {
  if (energy_map_) {
    if (energy_map_->is_queryable()) {
      return true;
    }
  }
  return false;
}

std::shared_ptr<VisitModelInner> VisitModelInner::deserialize(
    std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

void VisitModelInner::serialize_visit_model_inner_(std::ostream& ostr) const {
  feasst_serialize_version(8178, ostr);
  feasst_serialize(energy_, ostr);
  feasst_serialize(cutoff_index_, ostr);
  feasst_serialize(cutoff_outer_index_, ostr);
  feasst_serialize_fstdr(energy_map_, ostr);
}

VisitModelInner::VisitModelInner(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 8177 || version == 8178,
    "unrecognized verison: " << version);
  feasst_deserialize(&energy_, istr);
  feasst_deserialize(&cutoff_index_, istr);
  feasst_deserialize(&cutoff_outer_index_, istr);
  // HWH for unknown reasons, this function template does not work.
  { int existing;
    istr >> existing;
    if (existing != 0) {
      energy_map_ = energy_map_->deserialize(istr);
    }
  }
}

void VisitModelInner::query_ixn(const Select& select) {
  ASSERT(select.num_particles() <= 1, "not implemented");
  for (int pindex = 0; pindex < select.num_particles(); ++pindex) {
    const int part1_index = select.particle_index(pindex);
    for (const int site1_index : select.site_indices(pindex)) {
      energy_ += energy_map_->energy(part1_index, site1_index);
    }
  }
}

void VisitModelInner::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_visit_model_inner_(ostr);
}

void VisitModelInner::precompute(Configuration * config) {
  if (energy_map_) {
    energy_map_->precompute(config);
  }
  cutoff_index_ = config->model_params().index("cutoff");
  cutoff_outer_index_ = config->model_params().index("cutoff_outer");
  //TRACE("cutoff_outer_index " << cutoff_outer_index_);
}

}  // namespace feasst
