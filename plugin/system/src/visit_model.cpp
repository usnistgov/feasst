#include <vector>
#include "system/include/visit_model.h"
#include "system/include/model_two_body.h"
#include "system/include/model_one_body.h"
#include "system/include/select_list.h"
#include "math/include/constants.h"

namespace feasst {

void VisitModel::compute(
    const ModelOneBody& model,
    const ModelParams& model_params,
    Configuration * config,
    const int group_index) {
  zero_energy();
  const Select& selection = config->group_selects()[group_index];
  for (int select_index = 0;
       select_index < selection.num_particles();
       ++select_index) {
    const int part_index = selection.particle_index(select_index);
    const Particle& part = config->select_particle(part_index);
    for (int site_index : selection.site_indices(select_index)) {
      const Site& site = part.site(site_index);
      energy_ += model.energy(site, config, model_params);
    }
  }
}

void VisitModel::compute(
    const ModelOneBody& model,
    const ModelParams& model_params,
    const Select& selection,
    Configuration * config,
    const int group_index) {
  ASSERT(group_index == 0, "not implemented because redundant to selection");
  zero_energy();
  DEBUG("HWH: add wrapping of site positions");
  for (int select_index = 0; select_index < selection.num_particles(); ++select_index) {
    const int particle_index = selection.particle_index(select_index);
    const Particle& part = config->select_particle(particle_index);
    for (int site_index : selection.site_indices(select_index)) {
      const Site& site = part.site(site_index);
      energy_ += model.energy(site, config, model_params);
    }
  }
}

void VisitModel::compute(
    const ModelTwoBody& model,
    const ModelParams& model_params,
    Configuration * config,
    const int group_index) {
  zero_energy();
  const Domain& domain = config->domain();
  init_relative_(domain, &relative_);
  const Select& selection = config->group_selects()[group_index];
  for (int select1_index = 0;
       select1_index < selection.num_particles() - 1;
       ++select1_index) {
    const int part1_index = selection.particle_index(select1_index);
    for (int select2_index = select1_index + 1;
         select2_index < selection.num_particles();
         ++select2_index) {
      const int part2_index = selection.particle_index(select2_index);
      for (int site1_index : selection.site_indices(select1_index)) {
        for (int site2_index : selection.site_indices(select2_index)) {
          inner()->compute(part1_index, site1_index, part2_index, site2_index,
                           config, model_params, model, &relative_);
        }
      }
    }
  }
  set_energy(inner()->energy());
}

void VisitModel::compute(
    const ModelTwoBody& model,
    const ModelParams& model_params,
    const Select& selection,
    Configuration * config,
    const int group_index) {
  DEBUG("visiting model");
  zero_energy();
  const Domain& domain = config->domain();
  init_relative_(domain, &relative_);
  const Select& select_all = config->group_selects()[group_index];
  // HWH implement multi-particle selection by sorting group selection
  // for particles that are in both selectiona nd group_index.
  // treat those particles separately so no double counting.
  // then remove the part1 != part2 check
  ASSERT(selection.num_particles() == 1, "for multiparticle selections " <<
    "implement a separate loop for particles in both group and selection. " <<
    "Select: " << selection.str());
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
        for (int site1_index : selection.site_indices(select1_index)) {
          TRACE("site1_index " << site1_index);
          for (int site2_index : select_all.site_indices(select2_index)) {
            TRACE("index: " << part1_index << " " << part2_index << " " <<
                  site1_index << " " << site2_index);
            inner()->compute(part1_index, site1_index, part2_index, site2_index,
                             config, model_params, model, &relative_);
          }
        }
      }
    }
  }
  set_energy(inner()->energy());
}

void VisitModel::check_energy(
    const Model& model,
    Configuration * config,
    const int group_index) {
  model.compute(group_index, config, this);
  const double en_group = energy();

  // select each particle and compare half the sum with the whole
  SelectList select;
  double en_select = 0;
  const int num = config->num_particles(group_index);
  for (int part = 0; part < num; ++part) {
    select.particle(part, *config, group_index);
    model.compute(select, group_index, config, this);
    TRACE("part " << part << " en " << energy());
    en_select += 0.5*energy();
  }
  ASSERT(std::abs(en_group - en_select) < num*num*1e-15, "Error with " <<
    "visitor implementation. The energy of " <<
    MAX_PRECISION << "group(" << group_index << "): " << en_group << " "
    "is not consistent with half the sum of the energies of the selected " <<
    "particles: " << en_select << ". The difference is: " <<
    en_group - en_select << " with tolerance: " << num*num*1e-15);
}

void VisitModelInner::compute(
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
  if (squared_distance <= cutoff*cutoff) {
    energy_ += model.energy(squared_distance, type1, type2, model_params);
    TRACE("indices " << part1_index << " " << site1_index << " " <<
      part2_index << " " << site2_index);
    TRACE("energy " << energy_);
  }
}

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

class MapVisitModel {
 public:
  MapVisitModel() {
    VisitModel().deserialize_map()["VisitModel"] =
      std::make_shared<VisitModel>();
  }
};

static MapVisitModel mapper_visit_model_ = MapVisitModel();

std::map<std::string, std::shared_ptr<VisitModel> >& VisitModel::deserialize_map() {
  static std::map<std::string, std::shared_ptr<VisitModel> >* ans =
     new std::map<std::string, std::shared_ptr<VisitModel> >();
  return *ans;
}

void VisitModel::serialize_visit_model_(std::ostream& ostr) const {
  feasst_serialize_version(545, ostr);
  feasst_serialize(energy_, ostr);
  feasst_serialize_fstdr(inner_, ostr);
}

//void VisitModel::deserialize_visit_model_(std::istream& istr, std::shared_ptr<VisitModel> visitor) const {
VisitModel::VisitModel(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(545 == version, "mismatch: " << version);
  feasst_deserialize(&energy_, istr);
  // feasst_deserialize_fstdr(inner_, istr);
  { // for unknown reason, template function above does not work
    int existing;
    istr >> existing;
    if (existing != 0) {
      inner_ = inner_->deserialize(istr);
    }
  }
}

}  // namespace feasst
