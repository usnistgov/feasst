#include <vector>
#include "system/include/visit_model_inner.h"
#include "system/include/model_two_body.h"
#include "system/include/select_list.h"

namespace feasst {

class MapVisitModelInner {
 public:
  MapVisitModelInner() {
    VisitModelInner().deserialize_map()["VisitModelInner"] =
      std::make_shared<VisitModelInner>();
  }
};

static MapVisitModelInner mapper_visit_model_inner_ = MapVisitModelInner();

std::map<std::string, std::shared_ptr<VisitModelInner> >& VisitModelInner::deserialize_map() {
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
    const ModelTwoBody& model,
    const bool is_old_config,
    Position * relative,
    Position * pbc) {
  if (is_old_config && energy_map()) {
    DEBUG("using old map");
    query_ixn(part1_index, site1_index, part2_index, site2_index);
    return;
  }
  const Particle& part1 = config->select_particle(part1_index);
  const Site& site1 = part1.site(site1_index);
  clear_ixn(part1_index, site1_index, part2_index, site2_index);
  if (site1.is_physical()) {
    const Particle& part2 = config->select_particle(part2_index);
    const Site& site2 = part2.site(site2_index);
    if (site2.is_physical()) {
      config->domain().wrap_opt(site1.position(), site2.position(), relative, pbc, &squared_distance_);
      const int type1 = site1.type();
      const int type2 = site2.type();
      const double cutoff = model_params.mixed_cutoff()[type1][type2];
      if (squared_distance_ <= cutoff*cutoff) {
        const double energy = model.energy(squared_distance_, type1, type2, model_params);
        update_ixn(energy, part1_index, site1_index, part2_index, site2_index,
                   squared_distance_, pbc);
        //energy_ += model.energy(squared_distance_, type1, type2, model_params);
        TRACE("indices " << part1_index << " " << site1_index << " " <<
          part2_index << " " << site2_index);
        TRACE("energy " << energy_ << " " << energy);
      }
    }
  }
}

}  // namespace feasst
