#include <cmath>
#include <vector>
#include "utils/include/utils.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"
#include "configuration/include/select.h"
#include "configuration/include/configuration.h"
#include "configuration/include/domain.h"
#include "system/include/visit_model.h"
#include "system/include/model_two_body.h"
#include "system/include/model_one_body.h"

namespace feasst {

VisitModel::VisitModel(std::shared_ptr<VisitModelInner> inner) {
  set_inner(inner);
  energy_cutoff_ = -1;
}
VisitModel::VisitModel(argtype * args) {
  set_inner(VisitModelInner().factory(str("VisitModelInner", args, "VisitModelInner"), args));
  energy_cutoff_ = dble("energy_cutoff", args, -1);
  if (energy_cutoff_ != -1) {
    ASSERT(energy_cutoff_ > 1e10, "energy_cutoff:" << energy_cutoff_ <<
      " should be > 1e10 to avoid any trial with a chance of being accepted.");
  }
}
VisitModel::VisitModel(argtype args) : VisitModel(&args) {
  FEASST_CHECK_ALL_USED(args);
}
void VisitModel::compute(
    ModelOneBody * model,
    const ModelParams& model_params,
    Configuration * config,
    const int group_index) {
  zero_energy();
  const Domain& domain = config->domain();
  init_relative_(domain, &relative_, &pbc_);
  double r2;
  const Select& selection = config->group_selects()[group_index];
  for (int select_index = 0;
       select_index < selection.num_particles();
       ++select_index) {
    const int part_index = selection.particle_index(select_index);
    const Particle& part = config->select_particle(part_index);
    for (int site_index : selection.site_indices(select_index)) {
      const Site& site = part.site(site_index);
      if (site.is_physical()) {
        domain.wrap_opt(site.position(), origin_, &relative_, &pbc_, &r2);
        energy_ += model->energy(relative_, site, *config, model_params);
      }
    }
  }
}

void VisitModel::compute(
    ModelOneBody * model,
    const ModelParams& model_params,
    const Select& selection,
    Configuration * config,
    const int group_index) {
  ASSERT(group_index == 0, "not implemented because redundant to selection");
  zero_energy();
  const Domain& domain = config->domain();
  init_relative_(domain, &relative_, &pbc_);
  double r2;
  for (int sel_index = 0; sel_index < selection.num_particles(); ++sel_index) {
    const int particle_index = selection.particle_index(sel_index);
    const Particle& part = config->select_particle(particle_index);
    for (int site_index : selection.site_indices(sel_index)) {
      const Site& site = part.site(site_index);
      if (site.is_physical()) {
        domain.wrap_opt(site.position(), origin_, &relative_, &pbc_, &r2);
        energy_ += model->energy(relative_, site, *config, model_params);
      }
    }
  }
}

void VisitModel::compute(
    ModelTwoBody * model,
    const ModelParams& model_params,
    Configuration * config,
    const int group_index) {
  TRACE("VisitModel for TwoBody entire config");
  TRACE("VisitModelInner: " << get_inner_()->class_name());
  zero_energy();
  const Domain& domain = config->domain();
  init_relative_(domain, &relative_, &pbc_);
  TRACE("group index " << group_index);
  const Select& selection = config->group_selects()[group_index];
  TRACE("num p " << selection.num_particles());
  for (int select1_index = 0;
       select1_index < selection.num_particles() - 1;
       ++select1_index) {
    const int part1_index = selection.particle_index(select1_index);
    TRACE("part1_index " << part1_index);
    for (int select2_index = select1_index + 1;
         select2_index < selection.num_particles();
         ++select2_index) {
      const int part2_index = selection.particle_index(select2_index);
      for (int site1_index : selection.site_indices(select1_index)) {
        for (int site2_index : selection.site_indices(select2_index)) {
          get_inner_()->compute(part1_index, site1_index, part2_index,
            site2_index, config, model_params, model, false, &relative_, &pbc_);
          if ((energy_cutoff_ != -1) && (inner().energy() > energy_cutoff_)) {
            set_energy(inner().energy());
            return;
          }
        }
      }
    }
  }
  set_energy(inner().energy());
}

void VisitModel::compute(
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
  if (selection.num_particles() == 1) {
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
          for (const int site1_index : selection.site_indices(select1_index)) {
            TRACE("site1_index " << site1_index);
            for (const int site2_index : select_all.site_indices(select2_index)) {
              TRACE("index: " << part1_index << " " << part2_index << " " <<
                    site1_index << " " << site2_index);
              get_inner_()->compute(part1_index, site1_index,
                                    part2_index, site2_index,
                                    config, model_params, model,
                                    is_old_config,
                                    &relative_, &pbc_);
              if ((energy_cutoff_ != -1) && (inner().energy() > energy_cutoff_)) {
                set_energy(inner().energy());
                return;
              }
            }
          }
        }
      }
    }
  // If selection is more than one particle, skip those in selection
  // Calculate energy in two separate loops.
  } else {
    TRACE("more than one particle in selection");
    for (int select2_index = 0;
         select2_index < select_all.num_particles();
         ++select2_index) {
      const int part2_index = select_all.particle_index(select2_index);
      if (!find_in_list(part2_index, selection.particle_indices())) {
        for (int select1_index = 0;
             select1_index < selection.num_particles();
             ++select1_index) {
          const int part1_index = selection.particle_index(select1_index);
          TRACE("part1_index " << part1_index << " s " <<
                selection.particle_indices().size() << " " <<
                selection.site_indices().size());
          for (const int site1_index : selection.site_indices(select1_index)) {
            TRACE("site1_index " << site1_index);
            for (const int site2_index : select_all.site_indices(select2_index)) {
              TRACE("index: " << part1_index << " " << part2_index << " " <<
                    site1_index << " " << site2_index);
              get_inner_()->compute(part1_index, site1_index,
                                    part2_index, site2_index,
                                    config, model_params, model,
                                    is_old_config,
                                    &relative_, &pbc_);
              if ((energy_cutoff_ != -1) && (inner().energy() > energy_cutoff_)) {
                set_energy(inner().energy());
                return;
              }
            }
          }
        }
      }
    }

    // In the second loop, compute interactions between different particles in select.
    for (int select1_index = 0;
         select1_index < selection.num_particles() - 1;
         ++select1_index) {
      const int part1_index = selection.particle_index(select1_index);
      TRACE("sel1 " << select1_index << " part1_index " << part1_index << " s "
            << selection.particle_indices().size() << " " <<
            selection.site_indices().size());
      for (int select2_index = select1_index + 1;
           select2_index < selection.num_particles();
           ++select2_index) {
        const int part2_index = selection.particle_index(select2_index);
        if (part1_index != part2_index) {
          TRACE("sel2 " << select2_index << " part2_index " << part2_index);
          for (const int site1_index : selection.site_indices(select1_index)) {
            TRACE("site1_index " << site1_index);
            for (const int site2_index : selection.site_indices(select2_index)) {
              TRACE("index: " << part1_index << " " << part2_index << " " <<
                    site1_index << " " << site2_index);
              get_inner_()->compute(part1_index, site1_index,
                                    part2_index, site2_index,
                                    config, model_params, model,
                                    is_old_config,
                                    &relative_, &pbc_);
              if ((energy_cutoff_ != -1) && (inner().energy() > energy_cutoff_)) {
                set_energy(inner().energy());
                return;
              }
            }
          }
        }
      }
    }
  }
  set_energy(inner().energy());
}

void VisitModel::check_energy(
    Model * model,
    Configuration * config,
    const int group_index) {
  TRACE("checking energy");
  model->compute(group_index, config, this);
  const double en_group = energy();

  // select each particle and compare half the sum with the whole
  double en_select = 0;
  const int num = config->num_particles(group_index);
  for (int part = 0; part < num; ++part) {
    Select select;
    select.add_particle(config->select_particle(part), part);
    model->compute(select, group_index, config, this);
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

class MapVisitModel {
 public:
  MapVisitModel() {
    VisitModel().deserialize_map()["VisitModel"] =
      std::make_shared<VisitModel>();
  }
};

static MapVisitModel mapper_visit_model_ = MapVisitModel();

std::map<std::string, std::shared_ptr<VisitModel> >&
    VisitModel::deserialize_map() {
  static std::map<std::string, std::shared_ptr<VisitModel> >* ans =
     new std::map<std::string, std::shared_ptr<VisitModel> >();
  return *ans;
}

std::shared_ptr<VisitModel> VisitModel::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

std::shared_ptr<VisitModel> VisitModel::factory(const std::string name, argtype * args) {
  DEBUG("name: " << name << ", args: " << str(*args));
  return template_factory(deserialize_map(), name, args);
}

void VisitModel::serialize_visit_model_(std::ostream& ostr) const {
  feasst_serialize_version(545, ostr);
  feasst_serialize(energy_, ostr);
  feasst_serialize(epsilon_index_, ostr);
  feasst_serialize(sigma_index_, ostr);
  feasst_serialize(cutoff_index_, ostr);
  feasst_serialize(charge_index_, ostr);
  feasst_serialize(energy_cutoff_, ostr);
  feasst_serialize_fstdr(inner_, ostr);
  feasst_serialize_fstobj(data_, ostr);
  feasst_serialize_fstobj(manual_data_, ostr);
}

VisitModel::VisitModel(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(545 == version, "mismatch: " << version);
  feasst_deserialize(&energy_, istr);
  feasst_deserialize(&epsilon_index_, istr);
  feasst_deserialize(&sigma_index_, istr);
  feasst_deserialize(&cutoff_index_, istr);
  feasst_deserialize(&charge_index_, istr);
  feasst_deserialize(&energy_cutoff_, istr);
  // feasst_deserialize_fstdr(inner_, istr);
  { // for unknown reason, template function above does not work
    int existing;
    istr >> existing;
    if (existing != 0) {
      inner_ = inner_->deserialize(istr);
    }
  }
  feasst_deserialize_fstobj(&data_, istr);
  feasst_deserialize_fstobj(&manual_data_, istr);
}

void VisitModel::compute(
    ModelThreeBody * model,
    const ModelParams& model_params,
    Configuration * config,
    const int group_index) {
  FATAL("not implemented");
}

void VisitModel::init_relative_(const Domain& domain, Position * relative,
                                Position * pbc) {
  if (relative->dimension() != domain.dimension()) {
    relative->set_vector(domain.side_lengths().coord());
    pbc->set_vector(domain.side_lengths().coord());
    origin_ = Position(domain.dimension());
  }
}

void VisitModel::compute(
    ModelOneBody * model,
    Configuration * config,
    const int group_index) {
  const ModelParams& model_params = config->model_params();
  compute(model, model_params, config, group_index);
}

void VisitModel::compute(
    ModelOneBody * model,
    const Select& selection,
    Configuration * config,
    const int group_index) {
  const ModelParams& model_params = config->model_params();
  compute(model, model_params, selection, config, group_index);
}

void VisitModel::compute(
    ModelTwoBody * model,
    Configuration * config,
    const int group_index) {
  const ModelParams& model_params = config->model_params();
  compute(model, model_params, config, group_index);
}

void VisitModel::compute(
    ModelTwoBody * model,
    const Select& selection,
    Configuration * config,
    const int group_index) {
  const ModelParams& model_params = config->model_params();
  compute(model, model_params, selection, config, group_index);
}

void VisitModel::compute(
    ModelThreeBody * model,
    Configuration * config,
    const int group_index) {
  const ModelParams& model_params = config->model_params();
  compute(model, model_params, config, group_index);
}

void VisitModel::synchronize_(const VisitModel& visit,
    const Select& perturbed) {
  data_ = visit.data();
  inner_->synchronize_(visit.inner(), perturbed);
}

void VisitModel::precompute(Configuration * config) {
  inner_->precompute(config);
  epsilon_index_ = config->model_params().index("epsilon");
  sigma_index_ = config->model_params().index("sigma");
  cutoff_index_ = config->model_params().index("cutoff");
  charge_index_ = config->model_params().index("charge");
}

}  // namespace feasst
