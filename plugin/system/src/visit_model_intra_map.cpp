#include <vector>
#include "utils/include/serialize.h"
#include "utils/include/utils.h"
#include "configuration/include/configuration.h"
#include "system/include/visit_model_intra_map.h"
#include "system/include/model_two_body.h"

namespace feasst {

VisitModelIntraMap::VisitModelIntraMap(argtype args) {
  class_name_ = "VisitModelIntraMap";
  exclude_bonds_ = boolean("exclude_bonds", &args, false);
  exclude_angles_ = boolean("exclude_angles", &args, false);
}

void VisitModelIntraMap::precompute(Configuration * config) {
  VisitModel::precompute(config);

  // if map has not been sized yet, include all interactions except self.
  if (include_map_.size() == 0) {
    include_map_.resize(config->num_particle_types());
    for (int ptype = 0; ptype < config->num_particle_types(); ++ptype) {
      const int num_sites = config->particle_type(ptype).num_sites();
      include_map_[ptype].resize(num_sites);
      for (int site1 = 0; site1 < num_sites; ++site1) {
        include_map_[ptype][site1].resize(num_sites, 1);
        include_map_[ptype][site1][site1] = 0;
      }
    }
  }

  if (exclude_bonds_) {
    for (int ptype = 0; ptype < config->num_particle_types(); ++ptype) {
      for (const Bond& bond : config->particle_type(ptype).bonds()) {
        include_map_[ptype][bond.site(0)][bond.site(1)] = 0;
        include_map_[ptype][bond.site(1)][bond.site(0)] = 0;
      }
    }
  }

  if (exclude_angles_) {
    for (int ptype = 0; ptype < config->num_particle_types(); ++ptype) {
      for (const Angle& angle : config->particle_type(ptype).angles()) {
        include_map_[ptype][angle.site(0)][angle.site(2)] = 0;
        include_map_[ptype][angle.site(2)][angle.site(0)] = 0;
      }
    }
  }
}

void VisitModelIntraMap::compute(
    ModelTwoBody * model,
    const ModelParams& model_params,
    const Select& selection,
    Configuration * config,
    const int group_index) {
  TRACE("intra map particle energy_of_selection");
  ASSERT(group_index == 0,
    "need to implement site1 loop filtering particles by group");
  zero_energy();
  const Domain& domain = config->domain();
  init_relative_(domain, &relative_, &pbc_);
  for (int sp1index = 0;
       sp1index < static_cast<int>(selection.particle_indices().size());
       ++sp1index) {
    const int part1_index = selection.particle_index(sp1index);
    TRACE("part1_index: " << part1_index);
    const Particle& part1 = config->select_particle(part1_index);
    const int part_type = part1.type();

    // loop through all unique pairs in select
    const std::vector<int>& select_sites = selection.site_indices()[sp1index];
    for (int s1_index = 0;
         s1_index < static_cast<int>(select_sites.size()); ++s1_index) {
      const int site1_index = select_sites[s1_index];
      for (int s2_index = s1_index;
           s2_index < static_cast<int>(select_sites.size()); ++s2_index) {
        const int site2_index = select_sites[s2_index];
        if (include_map_[part_type][site1_index][site2_index] == 1) {
          TRACE("sites: " << site1_index << " " << site2_index);
          get_inner_()->compute(part1_index, site1_index, part1_index,
            site2_index, config, model_params, model, false, &relative_,
            &pbc_);
        }
      }
    }

    // if selection is not all sites in particle, then also loop through all
    // interactions between sites in particle not in selection
    if (static_cast<int>(select_sites.size()) < part1.num_sites()) {
      for (int site1_index = 0; site1_index < part1.num_sites(); ++site1_index) {
        if (!find_in_list(site1_index, select_sites)) {
          for (const int site2_index : select_sites) {
            if (include_map_[part_type][site1_index][site2_index] == 1) {
              TRACE("sites: " << site1_index << " " << site2_index);
              get_inner_()->compute(part1_index, site1_index, part1_index,
                site2_index, config, model_params, model, false, &relative_,
                &pbc_);
            }
          }
        }
      }
    }
  }
  set_energy(inner().energy());
}

class MapVisitModelIntraMap {
 public:
  MapVisitModelIntraMap() {
    auto obj = MakeVisitModelIntraMap();
    obj->deserialize_map()["VisitModelIntraMap"] = obj;
  }
};

static MapVisitModelIntraMap mapper_ = MapVisitModelIntraMap();

std::shared_ptr<VisitModel> VisitModelIntraMap::create(std::istream& istr) const {
  return std::make_shared<VisitModelIntraMap>(istr);
}

VisitModelIntraMap::VisitModelIntraMap(std::istream& istr) : VisitModel(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(2958 == version, version);
  feasst_deserialize(&exclude_bonds_, istr);
  feasst_deserialize(&exclude_angles_, istr);
  feasst_deserialize(&include_map_, istr);
}

void VisitModelIntraMap::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_visit_model_(ostr);
  feasst_serialize_version(2958, ostr);
  feasst_serialize(exclude_bonds_, ostr);
  feasst_serialize(exclude_angles_, ostr);
  feasst_serialize(include_map_, ostr);
}

void VisitModelIntraMap::compute(
    ModelTwoBody * model,
    const ModelParams& model_params,
    Configuration * config,
    const int group_index) {
  compute(model, model_params, config->selection_of_all(),
          config, group_index);
}

}  // namespace feasst
