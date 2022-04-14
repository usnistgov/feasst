#include <sstream>
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "system/include/visit_model_cutoff_outer.h"
#include "system/include/model_two_body.h"
#include "system/include/model_one_body.h"

namespace feasst {

VisitModelCutoffOuter::VisitModelCutoffOuter(argtype * args) : VisitModel() {
  class_name_ = "VisitModelCutoffOuter";
}
VisitModelCutoffOuter::VisitModelCutoffOuter(argtype args) : VisitModelCutoffOuter(&args) {
  check_all_used(args);
}
VisitModelCutoffOuter::VisitModelCutoffOuter(std::shared_ptr<VisitModelInner> inner,
  argtype args) : VisitModelCutoffOuter(args) {
  set_inner(inner);
}

void VisitModelCutoffOuter::compute(
    ModelTwoBody * model,
    const ModelParams& model_params,
    const Select& selection,
    Configuration * config,
    const int group_index) {
  DEBUG("visiting model");
  zero_energy();
  const Domain& domain = config->domain();
  init_relative_(domain, &relative_, &pbc_);
  const Select& select_all = config->group_selects()[group_index];
  bool is_old_config = false;
  if (selection.trial_state() == 0 ||
      selection.trial_state() == 2) {
    is_old_config = true;
  }

  // If possible, query energy map of old configuration instead of pair loop
  if (is_old_config) {
    if (selection.num_particles() == 1) {
      if (get_inner_()->is_energy_map_queryable()) {
        get_inner_()->query_ixn(selection);
        set_energy(inner().energy());
        return;
      }
    }
  }

  // If only one particle in selection, simply exclude part1==part2
  ASSERT(selection.num_particles(), "only implemented for single particles");
  for (int select1_index = 0;
       select1_index < selection.num_particles();
       ++select1_index) {
    const int part1_index = selection.particle_index(select1_index);
    TRACE("part1_index " << part1_index << " s " <<
          selection.particle_indices().size() << " " <<
          selection.site_indices().size());
    for (int select2_index = 0;
         select2_index < select_all.num_particles();
         ++select2_index) {
      const int part2_index = select_all.particle_index(select2_index);
      if (part1_index != part2_index) {
        get_inner_()->set_skip_particle(false);
        for (const int site1_index : selection.site_indices(select1_index)) {
          TRACE("site1_index " << site1_index);
          for (const int site2_index : select_all.site_indices(select2_index)) {
            TRACE("index: " << part1_index << " " << part2_index << " " <<
                  site1_index << " " << site2_index);
            if (!inner().skip_particle()) {
              get_inner_()->compute(part1_index, site1_index,
                                    part2_index, site2_index,
                                    config, model_params, model,
                                    is_old_config,
                                    &relative_, &pbc_);
            }
          }
        }
      }
    }
  }
  set_energy(inner().energy());
}

class MapVisitModelCutoffOuter {
 public:
  MapVisitModelCutoffOuter() {
    auto obj = MakeVisitModelCutoffOuter();
    obj->deserialize_map()["VisitModelCutoffOuter"] = obj;
  }
};

static MapVisitModelCutoffOuter mapper_ = MapVisitModelCutoffOuter();

VisitModelCutoffOuter::VisitModelCutoffOuter(std::istream& istr) : VisitModel(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(2081 == version, version);
}

void VisitModelCutoffOuter::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_visit_model_(ostr);
  feasst_serialize_version(2081, ostr);
}

class MapCutoffOuter {
 public:
  MapCutoffOuter() {
    auto obj = std::make_shared<CutoffOuter>();
    obj->deserialize_map()["cutoff_outer"] = obj;
  }
};

static MapCutoffOuter mapper_patch_angle_ = MapCutoffOuter();

void CutoffOuter::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_model_param_(ostr);
  feasst_serialize_version(2498, ostr);
}

CutoffOuter::CutoffOuter(std::istream& istr) : ModelParam(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2498, "mismatch version: " << version);
}

}  // namespace feasst
