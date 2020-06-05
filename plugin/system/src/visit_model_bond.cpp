#include <vector>
#include "configuration/include/configuration.h"
#include "system/include/visit_model_bond.h"
#include "system/include/model_two_body.h"
#include "utils/include/serialize.h"

namespace feasst {

void VisitModelBond::compute(
    const ModelTwoBody& model,
    const ModelParams& model_params,
    const Select& selection,
    Configuration * config,
    const int group_index) {
  TRACE("intra particle energy_of_selection");
  ASSERT(group_index == 0, "need to implement site1 loop filtering particles by group");
  zero_energy();
  const Domain& domain = config->domain();
  init_relative_(domain, &relative_, &pbc_);
  prep_for_revert(selection);
  for (int select_index = 0;
       select_index < selection.num_particles();
       ++select_index) {
    const int part_index = selection.particle_index(select_index);
    const Particle& part = config->select_particle(part_index);
    const int part_type = part.type();
    for (const Bond& bond : config->particle_type(part_type).bonds()) {
      const int site0 = bond.site_indices().front();
      const int site1 = bond.site_indices().back();
      TRACE("sites " << site0 << " " << site1);
      bool exclude = false;
      if (selection.old_bond()) {
        if ( (site0 == selection.site_indices()[0][0] &&
              site1 == selection.old_bond()->site_indices()[0][0]) ||
             (site1 == selection.site_indices()[0][0] &&
              site0 == selection.old_bond()->site_indices()[0][0]) ) {
          exclude = true;
        }
      }
      if (!exclude) {
        get_inner_()->compute(part_index, site0, part_index, site1,
          config, model_params, model, false, &relative_, &pbc_);
      }
    }
    for (const Angle& angle : config->particle_type(part_type).angles()) {
      const int site0 = angle.site_indices().front();
      const int site2 = angle.site_indices().back();
      TRACE("sites " << site0 << " " << site2);
      get_inner_()->compute(part_index, site0, part_index, site2,
        config, model_params, model, false, &relative_, &pbc_);
    }
    // force inclusion of new bond
    if (selection.new_bond()) {
      get_inner_()->compute(part_index, selection.site_index(0, 0),
                            part_index, selection.new_bond()->site_index(0, 0),
        config, model_params, model, false, &relative_, &pbc_);
    }
  }
  set_energy(inner().energy());
}

class MapVisitModelBond {
 public:
  MapVisitModelBond() {
    VisitModelBond().deserialize_map()["VisitModelBond"] =
      std::make_shared<VisitModelBond>();
  }
};

static MapVisitModelBond mapper_ = MapVisitModelBond();

std::shared_ptr<VisitModel> VisitModelBond::create(std::istream& istr) const {
  return std::make_shared<VisitModelBond>(istr);
}

VisitModelBond::VisitModelBond(std::istream& istr) : VisitModel(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(895 == version, version);
}

void VisitModelBond::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_visit_model_(ostr);
  feasst_serialize_version(895, ostr);
}

void VisitModelBond::compute(
    const ModelTwoBody& model,
    const ModelParams& model_params,
    Configuration * config,
    const int group_index) {
  compute(model, model_params, config->selection_of_all(),
          config, group_index);
}

}  // namespace feasst
