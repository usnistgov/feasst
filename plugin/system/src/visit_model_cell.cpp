#include "utils/include/arguments.h"
#include "utils/include/io.h"
#include "utils/include/utils.h"
#include "utils/include/serialize.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/group.h"
#include "configuration/include/domain.h"
#include "configuration/include/model_params.h"
#include "configuration/include/configuration.h"
#include "system/include/ideal_gas.h"
#include "system/include/model_two_body.h"
#include "system/include/model_three_body.h"
#include "system/include/cells.h"
#include "system/include/visit_model_inner.h"
#include "system/include/visit_model_cell.h"

namespace feasst {

FEASST_MAPPER(VisitModelCell,);

VisitModelCell::VisitModelCell(argtype * args) : VisitModel(args) {
  class_name_ = "VisitModelCell";
  cells_ = std::make_shared<Cells>();
  min_length_ = str("min_length", args, "max_cutoff");
  if (used("cell_group_index", *args)) {
    group_index_ = integer("cell_group_index", args);
    ASSERT(!used("cell_group", *args),
      "do not use both cell_group_index and cell_group at the same time.");
  } else {
    group_index_ = 0;
    group_ = str("cell_group", args, "");
  }
  ASSERT(group_index_ >= 0, "invalid group_index: " << group_index_);
}
VisitModelCell::VisitModelCell(argtype args) : VisitModelCell(&args) {
  feasst_check_all_used(args);
}
VisitModelCell::VisitModelCell(std::shared_ptr<VisitModelInner> inner,
  argtype args) : VisitModelCell(args) {
  set_inner(inner);
}

// This slow implementation which constructs a new position is const required
// for use in Check
int VisitModelCell::cell_id(const Domain& domain,
                            const Position& position) const {
  Position scaled = position;
  domain.cartesian2scaled_wrap(position, &scaled);
  return cells_->id(scaled.coord());
}

// HWH note if there are problems with scaled coordinates here, it probably
// means there is an issue with wrapping. As currently implemented, translations
// automatically wrap. So if you're doing a test without them you might run
// into this issue.
int VisitModelCell::cell_id_opt_(const Domain& domain,
                                 const Position& position) {
  init_relative_(domain);
  Position * scaled = relative_.get();
  domain.cartesian2scaled_wrap(position, scaled);
  const int cellid = cells_->id(scaled->coord());
  return cellid;
}

double VisitModelCell::min_len_(const Configuration& config) const {
  double min_length = -1;
  if (min_length_ == "max_sigma") {
    min_length = config.model_params().select("sigma").mixed_max();
  } else if (min_length_ == "max_cutoff") {
    min_length = config.model_params().select("cutoff").mixed_max();
  } else {
    min_length = str_to_double(min_length_);
  }
  // increase min_length to account for triclinic cells
  const double factor = config.domain().min_side_length() /
                        config.domain().inscribed_sphere_diameter();
  DEBUG("factor " << factor << " new min_length " << factor*min_length);
  return factor*min_length;
}

void VisitModelCell::precompute(Configuration * config) {
  DEBUG("precomputing");
  VisitModel::precompute(config);
  ASSERT(config->domain().side_lengths().size() > 0,
    "cannot define cells before domain sides");
  // ASSERT(!config->domain().is_tilted(), "implement triclinic");
  // This error check was moved to MonteCarlo::add(Potential)
  if (!group_.empty()) {
    group_index_ = config->group_index(group_);
  }
  if (cells_->type() == -1) {
    rebuild_(*config);
    config->increment_num_cell_lists();
    init_relative_(config->domain());
    position_tracker_(config->group_select(group_index_), config);
  }
  check(*config);
}

void VisitModelCell::rebuild_(const Configuration& config) {
  DEBUG("rebuilding");
  const double min_length = min_len_(config);
  Cells cells;
  cells.create(min_length, config.domain().side_lengths().coord());
  cells.set_group(group_index_);
  if (cells_->num_total() == 0) {
    // if first initialize of cells
    cells.set_type(config.num_cell_lists());
  } else {
    cells.set_type(cells_->type());
  }
  if (cells.num_total() > 0) {
    cells_ = std::make_shared<Cells>(cells);
  } else {
    FATAL("Requested cell list rejected: min_length:" << min_length <<
          " did not meet requirements when the minimum domain side length " <<
          "is " << config.domain().min_side_length());
  }
  DEBUG("num cells " << cells_->num_total());
  DEBUG("volume " << config.domain().volume());
}

void VisitModelCell::change_volume(const double delta_volume, const int dimension, Configuration * config) {
  // HWH optimize check if rebuild is necessary
  bool rebuild = true;
  if (rebuild) {
    rebuild_(*config);
    DEBUG("position updates after change volume rebuild");
    position_tracker_(config->group_select(group_index_), config);
  }
}

void VisitModelCell::compute(
    ModelTwoBody * model,
    const ModelParams& model_params,
    Configuration * config,
    const int group_index) {
  zero_energy();
  const Domain& domain = config->domain();
  ASSERT(group_index == group_index_, "not equivalent");
  VisitModelInner * inner = get_inner_();
  init_relative_(domain);

  /*
    Loop index nomenclature
    ends in 1 or 2 to represent the pair
    cell -> id of cell
    select -> selection inside each cell
    select_index -> index of selection
    part -> particle
    part_index -> index of particle
    site -> site
   */

  // loop through neighboring cells where cell1 < cell2 only
  for (int cell1 = 0; cell1 < cells_->num_total(); ++cell1) {
    const Select& select1 = cells_->particles()[cell1];
    for (int cell2 : cells_->neighbor()[cell1]) {
      if (cell1 < cell2) {
        const Select& select2 = cells_->particles()[cell2];
        for (int select1_index = 0;
             select1_index < select1.num_particles();
             ++select1_index) {
          const int part1_index = select1.particle_index(select1_index);
          for (int select2_index = 0;
               select2_index < select2.num_particles();
               ++select2_index) {
            const int part2_index = select2.particle_index(select2_index);
            if (part1_index != part2_index) {
              for (int site1_index : select1.site_indices(select1_index)) {
                for (int site2_index : select2.site_indices(select2_index)) {
                  inner->compute(part1_index, site1_index, part2_index,
                                        site2_index, config, model_params,
                                        model, false, relative_.get(), pbc_.get());
                  if ((energy_cutoff() != -1) && (inner->energy() > energy_cutoff())) {
                    set_energy(inner->energy());
                    return;
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // loop through the same cell only
  for (int cell1 = 0; cell1 < cells_->num_total(); ++cell1) {
    const Select& select = cells_->particles()[cell1];
    for (int select1_index = 0;
         select1_index < select.num_particles() - 1;
         ++select1_index) {
      const int part1_index = select.particle_index(select1_index);
      for (int select2_index = select1_index + 1;
           select2_index < select.num_particles();
           ++select2_index) {
        const int part2_index = select.particle_index(select2_index);
        if (part1_index != part2_index) {
          for (int site1_index : select.site_indices(select1_index)) {
            for (int site2_index : select.site_indices(select2_index)) {
              inner->compute(part1_index, site1_index, part2_index,
                                    site2_index, config, model_params, model,
                                    false, relative_.get(), pbc_.get());
              if ((energy_cutoff() != -1) && (inner->energy() > energy_cutoff())) {
                set_energy(inner->energy());
                return;
              }
            }
          }
        }
      }
    }
  }
  set_energy(inner->energy());
}

void VisitModelCell::compute(
    ModelTwoBody * model,
    const ModelParams& model_params,
    const Select& selection,
    Configuration * config,
    const int group_index) {
  DEBUG("VisitModelCell sel " << selection.str());
  zero_energy();
  const Domain& domain = config->domain();
  ASSERT(group_index == group_index_, "not equivalent");
  init_relative_(domain);
  VisitModelInner * inner = get_inner_();
  const bool is_old_config = is_old_config_(selection);
  if (is_queryable_(selection, is_old_config, inner)) {
    return;
  }

  // If only one particle in selection, simply exclude part1==part2
  DEBUG("num particles in selection " << selection.num_particles());
  if (selection.num_particles() == 1) {
    for (int select1_index = 0;
         select1_index < selection.num_particles();
         ++select1_index) {
      const int part1_index = selection.particle_index(select1_index);
      const Particle& part1 = config->select_particle(part1_index);
      for (int site1_index : selection.site_indices(select1_index)) {
        const Site& site1 = part1.site(site1_index);
        const int cell1_index = cell_id_opt_(domain, site1.position());
        for (int cell2_index : cells_->neighbor()[cell1_index]) {
          const Select& cell2_parts = cells_->particles()[cell2_index];
          for (int select2_index = 0;
               select2_index < cell2_parts.num_particles();
               ++select2_index) {
            const int part2_index = cell2_parts.particle_index(select2_index);
            if (part1_index != part2_index) {
              TRACE("indices " <<
                    feasst_str(cell2_parts.site_indices(select2_index)));
              for (int site2_index : cell2_parts.site_indices(select2_index)) {
                TRACE("index: " << part1_index << " " << part2_index << " " <<
                     site1_index << " " << site2_index);
                inner->compute(part1_index, site1_index, part2_index,
                                      site2_index, config, model_params, model,
                                      is_old_config, relative_.get(), pbc_.get());
                if ((energy_cutoff() != -1) && (inner->energy() > energy_cutoff())) {
                  set_energy(inner->energy());
                  return;
                }
              }
            }
          }
        }
      }
    }

  // If selection is more than one particle, skip those in selection
  // Calculate energy between particles in selection in two separate loops.
  } else {
    for (int select1_index = 0;
         select1_index < selection.num_particles();
         ++select1_index) {
      const int part1_index = selection.particle_index(select1_index);
      const Particle& part1 = config->select_particle(part1_index);
      for (int site1_index : selection.site_indices(select1_index)) {
        const Site& site1 = part1.site(site1_index);
        const int cell1_index = cell_id_opt_(domain, site1.position());
        for (int cell2_index : cells_->neighbor()[cell1_index]) {
          const Select& cell2_parts = cells_->particles()[cell2_index];
          for (int select2_index = 0;
               select2_index < cell2_parts.num_particles();
               ++select2_index) {
            const int part2_index = cell2_parts.particle_index(select2_index);
            if (!find_in_list(part2_index, selection.particle_indices())) {
              TRACE("indices " <<
                    feasst_str(cell2_parts.site_indices(select2_index)));
              for (int site2_index : cell2_parts.site_indices(select2_index)) {
                TRACE("index: " << part1_index << " " << part2_index << " " <<
                     site1_index << " " << site2_index);
                inner->compute(part1_index, site1_index, part2_index,
                                      site2_index, config, model_params, model,
                                      is_old_config, relative_.get(), pbc_.get());
                if ((energy_cutoff() != -1) && (inner->energy() > energy_cutoff())) {
                  set_energy(inner->energy());
                  return;
                }
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

void VisitModelCell::compute(
    ModelThreeBody * model,
    const ModelParams& model_params,
    Configuration * config,
    const int group_index) {
  DEBUG("VisitModelCell ThreeBody whole config");
  zero_energy();
  const Domain& domain = config->domain();
  ASSERT(group_index == group_index_, "not equivalent");
  VisitModelInner * inner = get_inner_();
  ModelTwoBody * two_body = model->get_two_body();
  init_relative_(domain);
  int num_pair = 0;

  // loop through neighboring cells where cell1 < cell2 only
  for (int cell1 = 0; cell1 < cells_->num_total(); ++cell1) {
    const Select& select1 = cells_->particles()[cell1];
    for (int cell2 : cells_->neighbor()[cell1]) {
      if (cell1 < cell2) {
        const Select& select2 = cells_->particles()[cell2];
        for (int select1_index = 0;
             select1_index < select1.num_particles();
             ++select1_index) {
          const int part1_index = select1.particle_index(select1_index);
          for (int select2_index = 0;
               select2_index < select2.num_particles();
               ++select2_index) {
            const int part2_index = select2.particle_index(select2_index);
            if (part1_index != part2_index) {
              for (int site1_index : select1.site_indices(select1_index)) {
                for (int site2_index : select2.site_indices(select2_index)) {
                  inner->compute(part1_index, site1_index, part2_index,
                                 site2_index, config, model_params,
                                 two_body, false, relative_.get(), pbc_.get());
                  record_pair_(part1_index, site1_index, part2_index, site2_index, *relative_, &num_pair, inner);
                  if ((energy_cutoff() != -1) && (inner->energy() > energy_cutoff())) {
                    set_energy(inner->energy());
                    return;
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // loop through the same cell only
  for (int cell1 = 0; cell1 < cells_->num_total(); ++cell1) {
    const Select& select = cells_->particles()[cell1];
    for (int select1_index = 0;
         select1_index < select.num_particles() - 1;
         ++select1_index) {
      const int part1_index = select.particle_index(select1_index);
      for (int select2_index = select1_index + 1;
           select2_index < select.num_particles();
           ++select2_index) {
        const int part2_index = select.particle_index(select2_index);
        if (part1_index != part2_index) {
          for (int site1_index : select.site_indices(select1_index)) {
            for (int site2_index : select.site_indices(select2_index)) {
              inner->compute(part1_index, site1_index, part2_index,
                             site2_index, config, model_params, two_body,
                             false, relative_.get(), pbc_.get());
              record_pair_(part1_index, site1_index, part2_index, site2_index, *relative_, &num_pair, inner);
              if ((energy_cutoff() != -1) && (inner->energy() > energy_cutoff())) {
                set_energy(inner->energy());
                return;
              }
            }
          }
        }
      }
    }
  }
  pair_pair_(num_pair, model, model_params, config, NULL);
  DEBUG("computed en: " << inner->energy());
  set_energy(inner->energy());
}

void VisitModelCell::compute(
    ModelThreeBody * model,
    const ModelParams& model_params,
    const Select& selection,
    Configuration * config,
    const int group_index) {
  DEBUG("VisitModelCell ThreeBody selection");
  IdealGas ideal_gas;
  VisitModelInner * inner = get_inner_();
  ModelTwoBody * two_body = model->get_two_body();
  zero_energy();
  const Domain& domain = config->domain();
  ASSERT(group_index == group_index_, "not equivalent");
  init_relative_(domain);
  const Select& select_all = config->group_select(group_index);
  const bool is_old_config = is_old_config_(selection);

  // HWH: To optimize query from EnergyMapNeighbor:
  // EnergyMapNeighbor stored squared distance, not vector separation
  // 1. "old" config can use EnergyMapNeighbor to build pairs_, but pair_pair_ remains the same
  // 2. "new" config can use EnergyMapNeighbor only for site3, and pair_pair_ remains
  // 3. pair_pair_ may be further optimized if EnergyMap stores vector separations.
  //    Then there would be less need to store pairs_ in VisitModel (no need in old config, and new config use new)
  //    pair_pair_ alternative must be implemented in EnergyMapNeighbor derived class

  // If only one particle in selection, simply exclude part1==part2
  DEBUG("num particles in selection " << selection.num_particles());
  if (selection.num_particles() == 1) {
    ASSERT(selection.num_sites() == 1, "pair_pair_ sel search assumes 1 site");
    int num_pair = 0;
    for (int select1_index = 0;
         select1_index < selection.num_particles();
         ++select1_index) {
      const int part1_index = selection.particle_index(select1_index);
      const Particle& part1 = config->select_particle(part1_index);
      for (int site1_index : selection.site_indices(select1_index)) {
        const Site& site1 = part1.site(site1_index);
        const int cell1_index = cell_id_opt_(domain, site1.position());
        for (int cell2_index : cells_->neighbor()[cell1_index]) {
          const Select& cell2_parts = cells_->particles()[cell2_index];
          for (int select2_index = 0;
               select2_index < cell2_parts.num_particles();
               ++select2_index) {
            const int part2_index = cell2_parts.particle_index(select2_index);
            if (part1_index != part2_index) {
              TRACE("indices " <<
                    feasst_str(cell2_parts.site_indices(select2_index)));
              for (int site2_index : cell2_parts.site_indices(select2_index)) {
                TRACE("index: " << part1_index << " " << part2_index << " " <<
                     site1_index << " " << site2_index);
                inner->compute(part1_index, site1_index, part2_index,
                               site2_index, config, model_params, two_body,
                               is_old_config, relative_.get(), pbc_.get());
                record_pair_(part1_index, site1_index, part2_index, site2_index, *relative_, &num_pair, inner);
                if ((energy_cutoff() != -1) && (inner->energy() > energy_cutoff())) {
                  set_energy(inner->energy());
                  return;
                }

                // now find all pairs of those paired with the pair
                if (inner->interacted()) {
                  for (int cell3_index : cells_->neighbor()[cell2_index]) {
                    const Select& cell3_parts = cells_->particles()[cell3_index];
                    for (int select3_index = 0;
                         select3_index < cell3_parts.num_particles();
                         ++select3_index) {
                      const int part3_index = cell3_parts.particle_index(select3_index);
                      DEBUG("part3_index " << part3_index);
                      if ((part3_index != part1_index) &&
                          (part3_index != part2_index)) {
                        for (const int site3_index : cell3_parts.site_indices(select3_index)) {
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
      }
    }
    pair_pair_(num_pair, model, model_params, config, &selection);
  } else if (selection.num_particles() == select_all.num_particles()) {
    compute(model, model_params, config, group_index);
  } else {
    FATAL("not implemented");
  }
  set_energy(inner->energy());
}

void VisitModelCell::position_tracker_(const Select& select,
    Configuration * config) {
  for (int spindex = 0; spindex < select.num_particles(); ++spindex) {
    DEBUG("spindex " << spindex << " of " << select.num_particles());
    const int particle_index = select.particle_index(spindex);
    for (const int site_index : select.site_indices(spindex)) {
      DEBUG("site_index " << site_index << " of " << static_cast<int>(select.site_indices(spindex).size()));
      ASSERT(site_index >= 0, "index error");
      DEBUG("update cells");
      DEBUG("group " << cells_->group());
      const int group_index = cells_->group();
      ASSERT(group_index >= 0, "error");
      DEBUG("config " << config);
      DEBUG("particle_index " << particle_index);
      DEBUG("num particles " << config->num_particles());
      const Particle& part = config->select_particle(particle_index);
      ASSERT(config, "error");
      DEBUG("group_index: " << group_index);
      DEBUG("group selects size: " << config->group_selects().size());
      const Select& sel = config->group_select(group_index);
      DEBUG("sel " << sel.str());
      DEBUG("is group empty: " << sel.is_group_empty());
      const Group& group = config->group_select(group_index).group();
      if (group.is_in(part, particle_index)) {
        const Site& site = part.site(site_index);
        if (group.is_in(site)) {
          const int cell_new = cell_id_opt_(config->domain(), site.position());
          if (!one_site_select_) {
            one_site_select_ = std::make_shared<Select>();
          }
          if (one_site_select_->num_particles() == 0) {
            one_site_select_->add_site(0, 0);
          }
          one_site_select_->set_particle(0, particle_index);
          one_site_select_->set_site(0, 0, site_index);
          ParticleFactory * particles = config->get_particles_();
          Site * sitep = particles->get_particle(particle_index)->get_site(site_index);
          if (cells_->type() < site.num_cells()) {
            DEBUG(cells_->type());
            const int cell_old = site.cell(cells_->type());
            DEBUG("index " << particle_index << " " << site_index);
            DEBUG("new cell " << cell_new << " old cell " << cell_old);
//            DEBUG("before new cell set: " <<
//              particles->particle(particle_index).site(
//              site_index).property("cell0"));
            sitep->set_cell(cells_->type(), cell_new);
            DEBUG(one_site_select_->str());
            DEBUG(cells_->num_total());
            cells_->update(*one_site_select_, cell_new, cell_old);
          } else {
            sitep->add_cell(cell_new);
            DEBUG("adding to cell list cllnw "
              << cell_new << " si " << site_index);
            cells_->add(*one_site_select_, cell_new);
          }
        }
      }
    }
  }
}

void VisitModelCell::finalize(const Select& select, Configuration * config) {
  VisitModel::finalize(select, config);
  if (select.trial_state() == 2) {
    // remove particles from cell
    for (const int particle_index : select.particle_indices()) {
      // note: somewhat derivative of position_tracker
      const int group_index = cells_->group();
      const Particle& part = config->select_particle(particle_index);
      const Group& group = config->group_select(group_index).group();
      if (group.is_in(part, particle_index)) {
        for (int site_index = 0; site_index < part.num_sites(); ++site_index) {
          const Site& site = part.site(site_index);
          if (group.is_in(site)) {
            if (cells_->type() < site.num_cells()) {
              const int cell_old = site.cell(cells_->type());
              Select select;
              select.add_site(particle_index, site_index);
              cells_->remove(select, cell_old);
            }
          }
        }
      }
    }
  } else {
    position_tracker_(select, config);
  }
}

void VisitModelCell::check(const Configuration& config) const {
  VisitModel::check(config);
  // for each site in config that has cells, check that the cell is correct.
  int num_sites_in_cell = 0;
  for (const int part_index : config.selection_of_all().particle_indices()) {
    for (const Site& site : config.select_particle(part_index).sites()) {
      if (site.num_cells() > cells_->type()) {
        ++num_sites_in_cell;
        const int old_cell = site.cell(cells_->type());
        const int cur_cell = cell_id(config.domain(), site.position());
        ASSERT(old_cell == cur_cell,
          "old_cell: " << old_cell << " != cur_cell: " << cur_cell);
      }
    }
  }
  ASSERT(num_sites_in_cell == cells_->num_sites(),
    "num sites with cells: " << num_sites_in_cell << " != " <<
    cells_->num_sites());
}

const Cells& VisitModelCell::cells() const { return *cells_; }

VisitModelCell::VisitModelCell(std::istream& istr) : VisitModel(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(755 == version, version);
  feasst_deserialize(&min_length_, istr);
  feasst_deserialize(&group_index_, istr);
  feasst_deserialize(&group_, istr);
//  feasst_deserialize_fstobj(&opt_origin_, istr);
//  feasst_deserialize_fstobj(&opt_rel_, istr);
//  feasst_deserialize_fstobj(&opt_pbc_, istr);
// HWH for unknown reasons, this function template does not work.
  //feasst_deserialize(cells_, istr);
  {
    int existing;
    istr >> existing;
    if (existing != 0) {
      cells_ = std::make_shared<Cells>(istr);
    }
  }
}

void VisitModelCell::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_visit_model_(ostr);
  feasst_serialize_version(755, ostr);
  feasst_serialize(min_length_, ostr);
  feasst_serialize(group_index_, ostr);
  feasst_serialize(group_, ostr);
//  feasst_serialize_fstobj(opt_origin_, ostr);
//  feasst_serialize_fstobj(opt_rel_, ostr);
//  feasst_serialize_fstobj(opt_pbc_, ostr);
  feasst_serialize(cells_, ostr);
  DEBUG("size: " << ostr.tellp());
}

}  // namespace feasst
