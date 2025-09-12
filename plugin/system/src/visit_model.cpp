#include <cmath>
#include <vector>
#include "utils/include/io.h"
#include "utils/include/arguments.h"
#include "utils/include/utils.h"
#include "utils/include/serialize_extra.h"
#include "math/include/constants.h"
#include "math/include/position.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/group.h"
#include "configuration/include/select.h"
#include "configuration/include/configuration.h"
#include "configuration/include/domain.h"
#include "configuration/include/model_params.h"
#include "system/include/model_one_body.h"
#include "system/include/model_two_body.h"
#include "system/include/model_three_body.h"
#include "system/include/ideal_gas.h"
#include "system/include/visit_model_inner.h"
#include "system/include/visit_model.h"

namespace feasst {

VisitModel::VisitModel() : VisitModel(std::make_shared<VisitModelInner>()) {}
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
  feasst_check_all_used(args);
}
VisitModel::~VisitModel() {}

void VisitModel::compute(
    ModelOneBody * model,
    const ModelParams& model_params,
    Configuration * config,
    const int group_index) {
  zero_energy();
  const Domain& domain = config->domain();
  init_relative_(domain);
  double r2;
  const Select& selection = config->group_select(group_index);
  for (int select_index = 0;
       select_index < selection.num_particles();
       ++select_index) {
    const int part_index = selection.particle_index(select_index);
    const Particle& part = config->select_particle(part_index);
    for (int site_index : selection.site_indices(select_index)) {
      const Site& site = part.site(site_index);
      if (site.is_physical()) {
        domain.wrap_opt(site.position(), *origin_, relative_.get(), pbc_.get(), &r2);
        energy_ += model->energy(*relative_, site, *config, model_params);
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
  zero_energy();
  const Domain& domain = config->domain();
  init_relative_(domain);
  double r2;
  if (group_index == 0) {
    for (int sel_index = 0; sel_index < selection.num_particles(); ++sel_index) {
      const int particle_index = selection.particle_index(sel_index);
      const Particle& part = config->select_particle(particle_index);
      for (int site_index : selection.site_indices(sel_index)) {
        const Site& site = part.site(site_index);
        if (site.is_physical()) {
          domain.wrap_opt(site.position(), *origin_, relative_.get(),
                          pbc_.get(), &r2);
          energy_ += model->energy(*relative_, site, *config, model_params);
        }
      }
    }
  } else {
    const Group& grp = config->group_select(group_index).group();
    for (int sel_index = 0; sel_index < selection.num_particles(); ++sel_index) {
      const int particle_index = selection.particle_index(sel_index);
      const Particle& part = config->select_particle(particle_index);
      if (grp.is_in(part, particle_index)) {
        for (int site_index : selection.site_indices(sel_index)) {
          const Site& site = part.site(site_index);
          if (grp.is_in(site)) {
            if (site.is_physical()) {
              domain.wrap_opt(site.position(), *origin_, relative_.get(),
                              pbc_.get(), &r2);
              energy_ += model->energy(*relative_, site, *config, model_params);
            }
          }
        }
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
  VisitModelInner * inner = get_inner_();
  TRACE("VisitModelInner: " << inner->class_name());
  zero_energy();
  const Domain& domain = config->domain();
  init_relative_(domain);
  TRACE("group index " << group_index);
  const Select& selection = config->group_select(group_index);
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
          inner->compute(part1_index, site1_index, part2_index,
            site2_index, config, model_params, model, false, relative_.get(), pbc_.get());
          if ((energy_cutoff_ != -1) && (inner->energy() > energy_cutoff_)) {
            set_energy(inner->energy());
            return;
          }
        }
      }
    }
  }
  TRACE("computed en: " << inner->energy());
  set_energy(inner->energy());
}

void VisitModel::compute(
    ModelTwoBody * model,
    const ModelParams& model_params,
    const Select& selection,
    Configuration * config,
    const int group_index) {
  TRACE("VistModel selection " << selection.str());
  VisitModelInner * inner = get_inner_();
  zero_energy();
  const Domain& domain = config->domain();
  init_relative_(domain);
  const Select& select_all = config->group_select(group_index);
  const bool is_old_config = is_old_config_(selection);
  if (is_queryable_(selection, is_old_config, inner)) {
    TRACE("queried");
    return;
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
              inner->compute(part1_index, site1_index,
                                    part2_index, site2_index,
                                    config, model_params, model,
                                    is_old_config,
                                    relative_.get(), pbc_.get());
              if ((energy_cutoff_ != -1) && (inner->energy() > energy_cutoff_)) {
                set_energy(inner->energy());
                return;
              }
            }
          }
        }
      }
    }
  } else if (selection.is_equal(config->selection_of_all())) {
  //} else if (selection.num_particles() == config->num_particles()) {
    DEBUG("computing entire select");
    compute_between_selection(model, model_params, selection,
      config, is_old_config, relative_.get(), pbc_.get());

  // If selection is more than one particle but not all particles, skip those in selection
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
              inner->compute(part1_index, site1_index,
                                    part2_index, site2_index,
                                    config, model_params, model,
                                    is_old_config,
                                    relative_.get(), pbc_.get());
              if ((energy_cutoff_ != -1) && (inner->energy() > energy_cutoff_)) {
                set_energy(inner->energy());
                return;
              }
            }
          }
        }
      }
    }

    // In the second loop, compute interactions between different particles in select.
    compute_between_selection(model, model_params, selection,
      config, is_old_config, relative_.get(), pbc_.get());
  }
  set_energy(inner->energy());
}

void VisitModel::compute_between_selection(
    ModelTwoBody * model,
    const ModelParams& model_params,
    const Select& selection,
    Configuration * config,
    const bool is_old_config,
    Position * relative,
    Position * pbc) {
  DEBUG("VisitModel computing between selection " << selection.str());
  VisitModelInner * inner = get_inner_();
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
            inner->compute(part1_index, site1_index,
                                  part2_index, site2_index,
                                  config, model_params, model,
                                  is_old_config,
                                  relative_.get(), pbc_.get());
            if ((energy_cutoff_ != -1) && (inner->energy() > energy_cutoff_)) {
              set_energy(inner->energy());
              return;
            }
          }
        }
      }
    }
  }
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

FEASST_MAPPER(VisitModel,);

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

void VisitModel::record_pair_(const int part1_index, const int site1_index, const int part2_index, const int site2_index, const Position& rel, int * num_pair, VisitModelInner * inner) {
  if (inner->interacted()) {
    if (*num_pair >= static_cast<int>(pairs_.size())) {
      pairs_.resize(pairs_.size() + 1);
      pairs_[*num_pair] = pair3body();
    }
    pair3body * pair = &(pairs_[*num_pair]);
    pair->part1 = part1_index;
    pair->site1 = site1_index;
    pair->part2 = part2_index;
    pair->site2 = site2_index;
    pair->rel = rel;
    (*num_pair) = *num_pair + 1;
  }
}

void VisitModel::compute(
    ModelThreeBody * model,
    const ModelParams& model_params,
    Configuration * config,
    const int group_index) {
  DEBUG("VisitModel for ThreeBody entire config");
  VisitModelInner * inner = get_inner_();
  ModelTwoBody * two_body = model->get_two_body();
  DEBUG("VisitModelInner: " << inner->class_name());
  zero_energy();
  const Domain& domain = config->domain();
  init_relative_(domain);
  const Select& selection = config->group_select(group_index);
  int num_pair = 0;
  for (int select1_index = 0;
       select1_index < selection.num_particles() - 1;
       ++select1_index) {
    const int part1_index = selection.particle_index(select1_index);
    DEBUG("part1_index " << part1_index);
    for (int site1_index : selection.site_indices(select1_index)) {
      DEBUG("site1_index " << site1_index);
      for (int select2_index = select1_index + 1;
           select2_index < selection.num_particles();
           ++select2_index) {
        const int part2_index = selection.particle_index(select2_index);
        DEBUG("part2_index " << part2_index);
        for (int site2_index : selection.site_indices(select2_index)) {
          DEBUG("site2_index " << site2_index);
          inner->compute(part1_index, site1_index, part2_index,
            site2_index, config, model_params, two_body, false, relative_.get(), pbc_.get());
          record_pair_(part1_index, site1_index, part2_index, site2_index, *relative_, &num_pair, inner);
        }
      }
    }
  }
  pair_pair_(num_pair, model, model_params, config, NULL);
  DEBUG("computed en: " << inner->energy());
  set_energy(inner->energy());
}

void VisitModel::pair_pair_(const int num_pair,
    ModelThreeBody * model,
    const ModelParams& model_params,
    Configuration * config,
    const Select * sel) {
  VisitModelInner * inner = get_inner_();
  Position flip;
  DEBUG("num_pair " << num_pair);
  int part1, part2, part3;
  int site1, site2, site3;
  // loop through pairs of pairs to compute three body interactions
  for (int pair1 = 0; pair1 < num_pair - 1; ++pair1) {
    const pair3body& p1 = pairs_[pair1];
    const int part1_index = p1.part1;
    const int site1_index = p1.site1;
    const int part2_index = p1.part2;
    const int site2_index = p1.site2;
    DEBUG("pair1 " << pair1 << " p2 " << part2_index << " s2 " << site2_index);
    const Position * r12 = &(p1.rel);
    for (int pair2 = pair1 + 1; pair2 < num_pair; ++pair2) {
      const pair3body& p2 = pairs_[pair2];
      const int part11_index = p2.part1;
      const int site11_index = p2.site1;
      const int part22_index = p2.part2;
      const int site22_index = p2.site2;
      //DEBUG("pair2 " << pair2 << "p3 " << part3_index << " s3 " << site3_index);
      const Position * r13 = &(p2.rel);
      part1 = -1;
      const Position * r12t = const_cast<const Position *>(r12);
      const Position * r13t = const_cast<const Position *>(r13);
      // the order matters for the angle
      if ((part11_index == part1_index) && (site11_index == site1_index)) {
        part1 = part1_index;
        site1 = site1_index;
        part2 = part2_index;
        site2 = site2_index;
        part3 = part22_index;
        site3 = site22_index;
      } else if ((part22_index == part1_index) && (site22_index == site1_index)) {
        part1 = part1_index;
        site1 = site1_index;
        part2 = part2_index;
        site2 = site2_index;
        part3 = part11_index;
        site3 = site11_index;
        flip = *r13;
        flip.multiply(-1);
        r13t = const_cast<const Position *>(&flip);
      } else if ((part11_index == part2_index) && (site11_index == site2_index)) {
        part1 = part2_index;
        site1 = site2_index;
        part2 = part1_index;
        site2 = site1_index;
        part3 = part22_index;
        site3 = site22_index;
        flip = *r12;
        flip.multiply(-1);
        r12t = const_cast<const Position *>(&flip);
      } else if ((part22_index == part2_index) && (site22_index == site2_index)) {
        part1 = part2_index;
        site1 = site2_index;
        part2 = part1_index;
        site2 = site1_index;
        part3 = part11_index;
        site3 = site11_index;
      }
      if (part1 != -1) {
        if (!sel) {
          inner->compute3body(part1, site1, part2,
            site2, part3, site3, *r12t, *r13t,
            config, model_params, model, false);
        } else {
          // If select provided, check that atleast one site is in selection
          // assumes one particle in selection
          if (part1 == sel->particle_index(0) ||
              part2 == sel->particle_index(0) ||
              part3 == sel->particle_index(0)) {
            inner->compute3body(part1, site1, part2,
              site2, part3, site3, *r12t, *r13t,
              config, model_params, model, false);
          }
        }
      }
    }
  }
}

bool VisitModel::find_in_pair3body(const int ipart, const int isite, const int num_pair) {
  if (num_pair == 0) return false;
  for (int i = 0; i < num_pair; ++i) {
    const pair3body& p = pairs_[i];
    if (p.part2 == ipart &&
        p.site2 == isite) {
      return true;
    }
  }
  return false;
}

/*
 * Looping through a selection to find all 3 body interactions requires
 * 1. finding all pairs with selection
 * 2. finding all pairs of those paired with selection
 * 3. computing 3 body interactions
 */
void VisitModel::compute(
    ModelThreeBody * model,
    const ModelParams& model_params,
    const Select& selection,
    Configuration * config,
    const int group_index) {
  DEBUG("VisitModel for ThreeBody of selection");
  IdealGas ideal_gas;
  VisitModelInner * inner = get_inner_();
  ModelTwoBody * two_body = model->get_two_body();
  zero_energy();
  const Domain& domain = config->domain();
  init_relative_(domain);
  const Select& select_all = config->group_select(group_index);
  if (selection.num_particles() == 1) {
    ASSERT(selection.num_sites() == 1, "pair_pair_ sel search assumes 1 site");
    int num_pair = 0;
    for (int select1_index = 0;
         select1_index < selection.num_particles();
         ++select1_index) {
      const int part1_index = selection.particle_index(select1_index);
      DEBUG("part1_index " << part1_index);
      for (int site1_index : selection.site_indices(select1_index)) {
        DEBUG("site1_index " << site1_index);
        for (int select2_index = 0;
             select2_index < select_all.num_particles();
             ++select2_index) {
          const int part2_index = select_all.particle_index(select2_index);
          DEBUG("part2_index " << part2_index);
          if (part1_index != part2_index) {
            for (const int site2_index : select_all.site_indices(select2_index)) {
              DEBUG("site2_index " << site2_index);
              inner->compute(part1_index, site1_index, part2_index,
                site2_index, config, model_params, two_body, false, relative_.get(), pbc_.get());
              record_pair_(part1_index, site1_index, part2_index, site2_index, *relative_, &num_pair, inner);
              // now find all pairs of those paired with the pair
              // unless 3 was already found as a neighbor of 1
              if (inner->interacted()) {
                for (int select3_index = 0;
                     select3_index < select_all.num_particles();
                     ++select3_index) {
                  const int part3_index = select_all.particle_index(select3_index);
                  DEBUG("part3_index " << part3_index);
                  if ((part3_index != part1_index) &&
                      (part3_index != part2_index)) {
                    for (const int site3_index : select_all.site_indices(select3_index)) {
                      DEBUG("site3_index " << site3_index);
                      // check if 3 was already found as a neighbor of 1
                      if (!find_in_pair3body(part3_index, site3_index, num_pair)) {
                        inner->compute(part3_index, site3_index, part2_index,
                          site2_index, config, model_params, &ideal_gas, false, relative_.get(), pbc_.get());
                        record_pair_(part3_index, site3_index, part2_index, site2_index, *relative_, &num_pair, inner);
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
    pair_pair_(num_pair, model, model_params, config, &selection);
  } else if (selection.num_particles() == select_all.num_particles()) {
    compute(model, model_params, config, group_index);
  } else {
    FATAL("not implemented");
  }
  DEBUG("computed en: " << inner->energy());
  set_energy(inner->energy());
}

void VisitModel::init_relative_(const Domain& domain) {
  if (!relative_) {
    relative_ = std::make_shared<Position>();
    pbc_ = std::make_shared<Position>();
    origin_ = std::make_shared<Position>();
  }
  if (relative_->dimension() != domain.dimension()) {
    relative_->set_vector(domain.side_lengths().coord());
    pbc_->set_vector(domain.side_lengths().coord());
    origin_ = std::make_shared<Position>(domain.dimension());
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

void VisitModel::compute(
    ModelThreeBody * model,
    const Select& selection,
    Configuration * config,
    const int group_index) {
  const ModelParams& model_params = config->model_params();
  compute(model, model_params, selection, config, group_index);
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

void VisitModel::set_inner(const std::shared_ptr<VisitModelInner> inner) {
  inner_ = inner;
}

const VisitModelInner& VisitModel::inner() const {
  return const_cast<VisitModelInner&>(*inner_);
}

void VisitModel::zero_energy() {
  energy_ = 0.;
  inner_->set_energy(0.);
}

void VisitModel::revert(const Select& select) { inner_->revert(select); }

void VisitModel::finalize(const Select& select, Configuration * config) {
  inner_->finalize(select);
}

void VisitModel::check(const Configuration& config) const {
  inner_->check(config);
}

VisitModelInner * VisitModel::get_inner_() const { return inner_.get(); }

void VisitModel::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_visit_model_(ostr);
}

bool VisitModel::is_old_config_(const Select& selection) const {
  bool is_old_config = false;
  if (selection.trial_state() == 0 ||
      selection.trial_state() == 2) {
    is_old_config = true;
  }
  DEBUG("is_old_config " << is_old_config);
  return is_old_config;
}

bool VisitModel::is_queryable_(const Select& selection, const bool is_old_config, VisitModelInner * inner) {
  if (is_old_config) {
    if (selection.num_particles() == 1) {
      if (inner->is_energy_map_queryable()) {
        inner->query_ixn(selection);
        TRACE("en: " << inner->energy());
        set_energy(inner->energy());
        return true;
      }
    }
  }
  return false;
}

}  // namespace feasst
