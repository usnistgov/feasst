#include <fstream>
#include <sstream>
#include "configuration/include/configuration.h"
#include "configuration/include/file_xyz.h"
#include "configuration/include/domain.h"
#include "utils/include/utils.h"
#include "math/include/utils_math.h"
#include "utils/include/debug.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"

namespace feasst {

Configuration::Configuration(const argtype& args) {
  domain_ = std::make_shared<Domain>();
  particle_types_.unique_particles();
  unique_types_.unique_types();
  // reset_unique_indices_();
  add(MakeGroup());  // add empty group which represents all particles
  args_.init(args);

  DEBUG("parse physical constants");
  if (args_.key("physical_constants").used()) {
    std::stringstream ss(args_.str());
    set_physical_constants(MakeCODATA2014()->deserialize(ss));
  }

  std::string start;

//  DEBUG("parse domain");
//  if (args_.key("cubic_box_length").used()) {
//    domain_ = MakeDomain({{"cubic_box_length", args_.str()}});
//  }

  DEBUG("parse types");
  // if only one particle type, allow drop the subscript
  start.assign("particle_type");
  if (args_.key(start).used()) {
    add_particle_type(args_.str());
  } else {
    int type = num_particle_types();
    std::stringstream key;
    key << start << type;
    while (args_.key(key.str()).used()) {
      add_particle_type(args_.str());
      ++type;
      ASSERT(type < 1e8, "type(" << type << ") is very high. Infinite loop?");
      key.str("");
      key << start << type;
    }
  }

  init_wrap(args_.key("wrap").dflt("true").boolean());
  args_.check_all_used();
}

void Configuration::add_particle_type(const std::string file_name,
    const std::string append) {
  DEBUG("adding type: " << file_name);
  ASSERT(particles_.num() == 0, "types cannot be added after particles");
  particle_types_.add(file_name);
  unique_types_.add(file_name);
  ghosts_.push_back(Select());
  ASSERT(ghosts_.back().is_group_empty(), "no ghosts in brand new type");
  ASSERT(!find_in_list(file_name + append, type_to_file_),
    "file_name(" << file_name << ") already provided.");
  type_to_file_.push_back(file_name + append);
  num_particles_of_type_.push_back(0);
}

void Configuration::add_(const Particle particle) {
  Particle part = particle;
  particles_.add(part);
  for (Select& select : group_selects_) {
    add_to_selection_(particles_.num() - 1, &select);
  }
  position_tracker_(particles_.num() - 1);
}

void Configuration::add_non_ghost_particle_of_type(const int type) {
  Particle part = particle_types_.particle(type);
  part.erase_bonds();
  add_(part);
  newest_particle_index_ = particles_.num() - 1;
  ++num_particles_of_type_[type];
}

void Configuration::add_particle_of_type(const int type) {
  DEBUG("adding type: " << type);
  ASSERT(type < num_particle_types(), "type(" << type << ") is not allowed "
    << "when there are only " << num_particle_types() << " particle types");
  // decide whether to add a ghost particle or a new one
  DEBUG("numg " << ghosts_[type].num_particles());
  if (ghosts_[type].num_particles() == 0) {
    add_non_ghost_particle_of_type(type);
  } else {
    const int index = ghosts_[type].particle_index(0);
    ghosts_[type].remove_particle(index);
    for (Select& select : group_selects_) {
      add_to_selection_(index, &select);
    }
    newest_particle_index_ = index;
    ++num_particles_of_type_[type];
  }
}

void Configuration::remove_particle_(const int particle_index) {
  // HWH optimization: this randomly generated string is expensive.
  // reset_unique_indices_();
  const int type = particles_.particle(particle_index).type();
  --num_particles_of_type_[type];
  ghosts_[type].add_particle(select_particle(particle_index), particle_index);
  DEBUG("type " << type);
  DEBUG("particle index " << particle_index);
  DEBUG("num particles " << num_particles());
  for (Select& select : group_selects_) {
    select.remove_particle(particle_index);
  }

  // remove particle_index from cell
  // note: somewhat derivative of position_tracker
  for (int cell_index = 0;
       cell_index < static_cast<int>(domain().cells().size());
       ++cell_index) {
    const Cells& cells = domain().cells()[cell_index];
    const int group_index = cells.group();
    const Particle& part = select_particle(particle_index);
    const Group& group = group_selects_[group_index].group();
    if (group.is_in(part)) {
      for (int site_index = 0; site_index < part.num_sites(); ++site_index) {
        const Site& site = part.site(site_index);
        if (group.is_in(site)) {
          if (cells.type() < site.num_cells()) {
            const int cell_old = site.cell(cells.type());
            Select select;
            select.add_site(particle_index, site_index);
            domain_->remove_from_cell_list(cell_index, select, cell_old);
          }
        }
      }
    }
  }
}

void Configuration::remove_particles(const Select& selection) {
  ASSERT(selection.num_particles() > 0, "no selection");
  for (int index = selection.num_particles()  - 1;
       index >= 0;
       --index) {
    remove_particle_(selection.particle_index(index));
  }
}

void Configuration::remove_particle(const Select& selection) {
  ASSERT(selection.num_particles() == 1, "assumes 1 particle selected");
  remove_particles(selection);
}

void Configuration::displace_particle_(const int particle_index,
                                       const Position &displacement) {
  particles_.displace(particle_index, displacement);
  position_tracker_(particle_index);
}

void Configuration::displace_site_(const int particle_index,
                                   const int site_index,
                                   const Position &displacement) {
  Position pos = displacement;
  pos.add(select_particle(particle_index).site(site_index).position());
  particles_.replace_position(particle_index, site_index, pos);
  position_tracker_(particle_index, site_index);
}

void Configuration::displace_particles(const Select& selection,
                                       const Position &displacement) {
  ASSERT(selection.num_particles() > 0, "no selection");
  for (int particle_index : selection.particle_indices()) {
    displace_particle_(particle_index, displacement);
  }
}

void Configuration::displace_particle(const Select& selection,
                                      const Position &displacement) {
  ASSERT(selection.num_particles() == 1, "assumes 1 particle selected");
  displace_particles(selection, displacement);
}

void Configuration::replace_position(const Select& select,
                                     const Particle& replacement) {
  ASSERT(select.num_particles() == 1, "only one particle can be selected");
  replace_position_(select.particle_index(0), replacement);
}

void Configuration::replace_position_(const int particle_index,
                                      const Particle& replacement) {
  ASSERT(select_particle(particle_index).num_sites() ==
         replacement.num_sites(), "size error");
  particles_.replace_position(particle_index, replacement);
  position_tracker_(particle_index);
}

void Configuration::replace_position_(const int particle_index,
                                      const int site_index,
                                      const Position& replacement) {
  particles_.replace_position(particle_index, site_index, replacement);
  position_tracker_(particle_index, site_index);
}

void Configuration::replace_position_(const int particle_index,
                                      const Position& replacement) {
  particles_.replace_position(particle_index, replacement);
  /// HWH no position_tracker_ for just particle positions.
}

void Configuration::add(std::shared_ptr<Group> group, std::string name) {
  ASSERT(group->is_empty() || particle_types_.num() != 0,
    "add groups after particle types");
  if (name == "") {
    std::stringstream ss;
    ss << num_groups();
    name.assign(ss.str());
  }
  Select group_select;
  group->add_property(name, 0.);
  group_select.set_group(group);
  init_selection_(&group_select);
  group_selects_.push_back(group_select);
}

void Configuration::position_tracker_(const int particle_index,
                                      const int site_index) {
  ASSERT(site_index >= 0, "index error");
  DEBUG("update cells");
  if (domain_) {
    for (int cell_index = 0;
         cell_index < static_cast<int>(domain().cells().size());
         ++cell_index) {
      DEBUG("cell index " << cell_index);
      const Cells& cells = domain().cells()[cell_index];
      DEBUG("group " << cells.group());
      const int group_index = cells.group();
      ASSERT(group_index >= 0, "error");
      const Particle& part = select_particle(particle_index);
      const Group& group = group_selects_[group_index].group();
      if (group.is_in(part)) {
        const Site& site = part.site(site_index);
        if (group.is_in(site)) {
          const int cell_new = domain_->cell_id(site.position(), cells);
          if (one_site_select_.num_particles() == 0) {
            one_site_select_.add_site(0, 0);
          }
          one_site_select_.set_particle(0, particle_index);
          one_site_select_.set_site(0, 0, site_index);
          Site * sitep = particles_.get_particle(particle_index)->get_site(site_index);
          if (cells.type() < site.num_cells()) {
            const int cell_old = site.cell(cells.type());
            DEBUG("index " << particle_index << " " << site_index);
            DEBUG("new cell " << cell_new << " old cell " << cell_old);
            DEBUG("before new cell set: " <<
              particles_.particle(particle_index).site(
              site_index).property("cell0"));
            sitep->set_cell(cells.type(), cell_new);
            domain_->update_cell_list(cell_index, one_site_select_, cell_new,
                                      cell_old);
          } else {
            sitep->add_cell(cell_new);
            DEBUG("adding to cell list " << cell_index << " cllnw "
              << cell_new << " si " << site_index);
            domain_->add_to_cell_list(cell_index, one_site_select_, cell_new);
          }
        }
      }
    }
  }

  DEBUG("update selection");
  for (const Select& select : group_selects_) {
    ASSERT(!select.group().is_spatial(), "implement updating of groups");
  }
  /// HWH update neighbors?
}

void Configuration::position_tracker_(const int particle_index) {
  for (int site_index = 0;
       site_index < select_particle(particle_index).num_sites();
       ++site_index) {
    position_tracker_(particle_index, site_index);
  }
}

void Configuration::position_tracker_() {
  for (int index : selection_of_all().particle_indices()) {
    position_tracker_(index);
  }
}

void Configuration::set(std::shared_ptr<Domain> domain) {
  domain_ = domain;
  position_tracker_();
}

/// HWH add check .. domain, positions, particles, etc
void Configuration::check() const {
  particles_.check();
  particle_types_.check();
  unique_types_.check();
  selection_of_all().check();

  ASSERT(particle_types_.num() == num_particle_types(), "er");
  ASSERT(unique_types_.num_sites() == num_site_types(), "er");

  // check that the first group is all particles in the configuration.
  ASSERT(static_cast<int>(group_selects_[0].num_particles()) == num_particles(),
    "The number of particles in the first group(" <<
    group_selects_[0].num_particles() << ") is not equal to the number of " <<
    "particles: " << num_particles());
  for (int index = 0; index < num_particles(); ++index) {
    ASSERT(static_cast<int>(group_selects_[0].site_indices(index).size()) ==
      particle(index).num_sites(), "size error");
  }

  if (domain_) {
    // check that number of particles in cell list is number in selection.
    for (const Cells& cells : domain().cells()) {
      const int group_index = cells.group();
      const Select& select = group_selects_[group_index];
      ASSERT(select.num_sites() == cells.num_sites(),
        "sites in group(" << select.num_sites() << ") is not equal to sites " <<
        "in cell(" << cells.num_sites() << ") for group: " << group_index);
    }
  }

  // check number of particle types
  ASSERT(
    static_cast<int>(num_particles_of_type_.size()) == num_particle_types(),
    "size error");

  // check that a particle is not simultaneously a ghost and a real particle
  for (const Select& ghost : ghosts_) {
    ASSERT(!group_selects_[0].is_overlap(ghost),
      "ghost particle cannot also be real");
  }

  // check that a particle is not listed as a ghost twice
  for (const Select ghost : ghosts_) {
    ASSERT(!has_duplicate(ghost.particle_indices()),
      "the same particle cannot be listed as a ghost twice");
  }

  model_params().check();
}

bool Configuration::are_all_sites_physical() const {
  for (const Particle& part : particles_.particles()) {
    for (const Site& site : part.sites()) {
      if (!site.is_physical()) {
        return false;
      }
    }
  }
  return true;
}

// void Configuration::check_id_(const Select& select) const {
//   check_id_(select.unique_id());
// }
// void Configuration::check_id_(const std::string id) const {
//   ASSERT(id.empty() or id == unique_indices_,
//     "If selection is obtained from a configuration, it will have a unique" <<
//     "identifier attached to it. The id of the selection and the " <<
//     "configuration were fount to be inconsistent, likely due to an " <<
//     "expired selection (e.g., particles were removed after selection)");
// }

void Configuration::update_positions(
    const std::vector<std::vector<double> > coords) {
  ASSERT(static_cast<int>(coords.size()) == num_sites(), "size error");
  DEBUG("dimension: " << dimension());
  ASSERT(static_cast<int>(coords[0].size()) == dimension(), "size error");
  Position position;
  int iter_site = 0;
  for (int part_index : group_selects_[0].particle_indices()) {
    int num_site = 0;
    Particle part = select_particle(part_index);
    for (int site_index = 0;
         site_index < part.num_sites();
         ++site_index) {
      position.set_vector(coords[iter_site]);
      Site site = part.site(site_index);
      site.set_position(position);
      if (num_site == 0) {
        part.set_position(position);
      }
      part.set_site(site_index, site);
      ++num_site;
      ++iter_site;
    }
    replace_position_(part_index, part);
  }
}

void Configuration::displace(const Select& selection,
                             const Position &displacement) {
  for (const int& particle_index : selection.particle_indices()) {
    for (const int& site_index : selection.site_indices(particle_index)) {
      displace_site_(particle_index, site_index, displacement);
    }
  }
}

void Configuration::add_to_selection_(const int particle_index,
                                      Select * select) const {
  const Particle& part = select_particle(particle_index);
  const Group& group = select->group();
  if (group.is_in(part)) {
    select->add_particle(particle_index, group.site_indices(part));
  }
}

void Configuration::update_selection_(const int particle_index,
                                      Select * select) const {
  const Particle& part = select_particle(particle_index);
  const Group& group = select->group();
  if (group.is_in(part)) {
    select->add_particle(particle_index, group.site_indices(part), true);
  } else {
    select->remove_particle(particle_index);
  }
}

void Configuration::init_selection_(Select * group_select) const {
  if (num_particles() > 0) {
    for (int part_index : selection_of_all().particle_indices()) {
      add_to_selection_(part_index, group_select);
    }
  }
}

void Configuration::update_positions(const Select& select,
                                     const bool no_wrap,
                                     const bool no_exclude) {
  int pindex = 0;
  for (int particle_index : select.particle_indices()) {
    replace_position_(particle_index, select.particle_positions()[pindex]);
    int sindex = 0;
    for (int site_index : select.site_indices(pindex)) {
      // DEBUG(select.site_properties()[pindex][sindex].str());
      const std::vector<std::string> * excluded = &excluded_properties_;
      if (no_exclude) excluded = &excl_prop_non_usr_;
      replace_properties_(particle_index,
                          site_index,
                          select.site_properties()[pindex][sindex],
                          *excluded);
      replace_position_(particle_index,
                        site_index,
                        select.site_positions()[pindex][sindex]);
      ++sindex;
    }

    // check if particle needs to be wrapped, for whole-particle updates only.
    if (!no_wrap) {
      if (select.num_sites(pindex) ==
          select_particle(particle_index).num_sites()) {
        DEBUG("wrapping");
        wrap_particle(particle_index);
      }
    }
    ++pindex;
  }
}

int Configuration::particle_type_to_group(const int particle_type) const {
  int index;
  if (num_particle_types() == 1) {
    return 0;
  }
  DEBUG("group_stored.. " << feasst_str(group_store_particle_type_));
  DEBUG("num_groups " << num_groups());
  DEBUG("particle_type " << particle_type);
  if (!find_in_list(particle_type, group_store_particle_type_, &index)) {
    return -1;
  }
  return group_store_group_index_[index];
}

int Configuration::particle_type_to_group_create(const int particle_type) {
  const int grp = particle_type_to_group(particle_type);
  if (grp != -1) {
    return grp;
  }
  add(MakeGroup({{"add_particle_type", str(particle_type)}}));
  group_store_particle_type_.push_back(particle_type);
  group_store_group_index_.push_back(num_groups() - 1);
  const int index = static_cast<int>(group_store_group_index_.size()) - 1;
  DEBUG("index " << index);
  DEBUG("group_stored p.. " << feasst_str(group_store_particle_type_));
  DEBUG("group_stored g.. " << feasst_str(group_store_group_index_));
  return group_store_group_index_[index];
}

int Configuration::num_particles(const int group) const {
  TRACE("pn " << particles_.num());
  TRACE("gn " << num_ghosts_());
  TRACE("group " << group);
  if (group == 0) {
    TRACE("here " << particles_.num() - num_ghosts_());
    const int num = particles_.num() - num_ghosts_();
    ASSERT(num >= 0, "error");
    return num;
  }
  TRACE("here " << group_selects_[group].num_particles());
  const int num = group_selects_[group].num_particles();
  ASSERT(num >= 0, "error");
  return num;
}

int Configuration::num_sites(const int group) const {
  if (group != 0) FATAL("not implemented");

  int num_ghost_sites = 0;
  for (const Select& ghost : ghosts_) {
    num_ghost_sites += ghost.num_sites();
  }
  DEBUG("num ghost sites " << num_ghost_sites);
  return particles_.num_sites() - num_ghost_sites;
}

int Configuration::num_ghosts_() const {
  int num = 0;
  for (const Select& select : ghosts_) {
    num += select.num_particles();
  }
  return num;
}

void Configuration::revive(const Select& selection) {
  for (int particle_index : selection.particle_indices()) {
    DEBUG("reviving particle_index: " << particle_index);
    const Particle& part = select_particle(particle_index);
    const int type = part.type();
    ++num_particles_of_type_[type];
    DEBUG("ghost particles " << ghosts_[type].num_particles());
    ghosts_[type].remove_particle(particle_index);
    for (Select& select : group_selects_) {
      add_to_selection_(particle_index, &select);
    }
    position_tracker_(particle_index);
  }
}

const Particle Configuration::particle(const int index,
                                       const int group) const {
  const Select& select_group = group_selects_[group];
  const int particle_index = select_group.particle_index(index);
  Particle part = particles_.particle(particle_index);
  select_group.group().remove_sites(&part);
  return part;
}

void Configuration::add(std::shared_ptr<ModelParam> param) {
  for (const Particle& particle : unique_types().particles()) {
    param->add(particle);
  }
  unique_types_.add(param);
}

int Configuration::num_particles_of_type(const int type) const {
  ASSERT(type < num_particle_types(), "type: " << type << " >= num types: "
    << num_particle_types());
  if (type == -1) {
    return num_particles();
  }
  const int num = num_particles_of_type_[type];
  ASSERT(num >= 0, "unphysical number(" << num <<")");
  ASSERT(num <= num_particles(), "unphysical number of particles of type(" <<
    num << ") when number of particles is: " << num_particles());
  return num;
}

void Configuration::wrap_particle(const int particle_index) {
  if (wrap_) {
    const Position& part_position = select_particle(particle_index).position();
    const Position& pbc_shift = domain_->shift_opt(part_position);
    DEBUG("part_position " << part_position.str());
    DEBUG("pbc " << pbc_shift.str());
    if (pbc_shift.squared_distance() > NEAR_ZERO) {
      displace_particle_(particle_index, pbc_shift);
      DEBUG("new pos " << select_particle(particle_index).position().str());
    }
  }
}

void Configuration::set_selection_physical(const Select& select,
    const bool phys) {
  for (int sp_index = 0;
       sp_index < static_cast<int>(select.particle_indices().size());
       ++sp_index) {
    const std::vector<int>& site_indices = select.site_indices()[sp_index];
    for (int ss_index = 0;
         ss_index < static_cast<int>(site_indices.size());
         ++ss_index) {
      particles_.set_site_physical(
        select.particle_indices()[sp_index],
        site_indices[ss_index],
        phys);
    }
  }
}

bool Configuration::is_equal(const Configuration& configuration,
                             const double tolerance) const {
  // check particles/sites of non-ghosts.
  if (!selection_of_all().is_equal(configuration.selection_of_all())) {
    DEBUG("unequal selection");
    return false;
  }

  // check position of first particle.
  if (num_particles() > 0) {
    for (int pindex = 0; pindex < num_particles(); ++pindex) {
      const Particle p1 = particle(pindex);
      const Particle p2 = configuration.particle(pindex);
      if (!p1.position().is_equal(p2.position(), tolerance)) {
        DEBUG("pindex " << pindex);
        DEBUG("unequal positions: " << p1.position().str() << " vs "
          << p2.position().str());
        return false;
      }
      for (int is = 0; is < p1.num_sites(); ++is) {
        if (!p1.site(is).position().is_equal(p2.site(is).position(),
                                             tolerance)) {
          DEBUG("unequal site" << is << " positions: "
            << p1.site(is).position().str() << " vs "
            << p2.site(is).position().str());
          return false;
        }
      }
    }
  }
  return true;
}

void Configuration::add_excluded_property(const std::string name) {
  if (!find_in_list(name, excluded_properties_)) {
    excluded_properties_.push_back(name);
  }
}

int Configuration::max_sites_in_any_particle() const {
  int mx = 0;
  for (int type = 0; type < num_particle_types(); ++type) {
    mx = std::max(mx, particle_type(type).num_sites());
  }
  return mx;
}

void Configuration::set_site_type(const int particle_type,
                                  const int site,
                                  const int site_type) {
  for (int cell_index = 0;
       cell_index < static_cast<int>(domain().cells().size());
       ++cell_index) {
    const Cells& cells = domain().cells()[cell_index];
    const int group_index = cells.group();
    ASSERT(group_index == 0,
      "check if cell list needs to be updated with changing type");
  }
  for (const Select& group : group_selects_) {
    if (find_in_list(site_type, group.group().site_types())) {
      ERROR("check if groups need to be updated with changing type");
    }
  }
  particle_types_.set_site_type(particle_type, site, site_type);
  for (int particle = 0; particle < particles_.num(); ++particle) {
    if (particles_.particle(particle).type() == particle_type) {
      particles_.set_site_type(particle, site, site_type);
    }
  }
}

void Configuration::recenter_particle_positions(const int group_index) {
  const Select& group = group_selects_[group_index];
  for (int spindex = 0; spindex < group.num_particles(); ++spindex) {
    const int pindex = group.particle_index(spindex);
    Select sel(pindex, select_particle(pindex));
    sel.set_particle_position(0, sel.geometric_center());
    update_positions(sel);
  }
}

std::vector<int> Configuration::num_sites_of_type(
    const Select& selection) const {
  std::vector<int> count(num_site_types());
  for (int select_index = 0;
       select_index < selection.num_particles();
       ++select_index) {
    const int part_index = selection.particle_index(select_index);
    const Particle& part = select_particle(part_index);
    for (int site_index : selection.site_indices(select_index)) {
      const Site& site = part.site(site_index);
      if (site.is_physical()) {
        ++count[site.type()];
      }
    }
  }
  return count;
}

void Configuration::set_side_lengths(const Position& sides) {
  if (domain_) {
    domain_->set_side_lengths(sides);
  } else {
    WARN("site lengths were attempted to be set without a domain");
  }
}

std::string Configuration::status_header() const {
  std::stringstream ss;
  ss << domain().status_header();
  for (int type = 0; type < num_particle_types(); ++type) {
    ss << ",p" << type;
  }
  return ss.str();
}

std::string Configuration::status() const {
  std::stringstream ss;
  ss << domain().status();
  for (int type = 0; type < num_particle_types(); ++type) {
    ss << "," << num_particles_of_type(type);
  }
  return ss.str();
}

void Configuration::serialize(std::ostream& ostr) const {
  feasst_serialize_version(7199, ostr);
  particle_types_.serialize(ostr);
  unique_types_.serialize(ostr);
  particles_.serialize(ostr);
  feasst_serialize(domain_, ostr);
  feasst_serialize_fstobj(group_selects_, ostr);
  feasst_serialize(group_store_particle_type_, ostr);
  feasst_serialize(group_store_group_index_, ostr);
  feasst_serialize_fstobj(ghosts_, ostr);
  feasst_serialize(type_to_file_, ostr);
  feasst_serialize(num_particles_of_type_, ostr);
  feasst_serialize(wrap_, ostr);
  feasst_serialize_endcap("Configuration", ostr);
}

Configuration::Configuration(std::istream& istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 7199, "unrecognized version: " << version);
  particle_types_ = ParticleFactory(istr);
  unique_types_ = ParticleFactory(istr);
  particles_ = ParticleFactory(istr);
  // HWH for unknown reasons, this doesn't work
  // feasst_deserialize(domain_, istr);
  { int existing;
    istr >> existing;
    if (existing != 0) domain_ = std::make_shared<Domain>(istr);
  }
  feasst_deserialize_fstobj(&group_selects_, istr);
  feasst_deserialize(&group_store_particle_type_, istr);
  feasst_deserialize(&group_store_group_index_, istr);
  feasst_deserialize_fstobj(&ghosts_, istr);
  feasst_deserialize(&type_to_file_, istr);
  feasst_deserialize(&num_particles_of_type_, istr);
  feasst_deserialize(&wrap_, istr);
  feasst_deserialize_endcap("Configuration", istr);
}

void Configuration::copy_particles(const Configuration& config,
    const bool add_missing) {
  for (int ipart = 0; ipart < num_particles(); ++ipart) {
    const int type1 = particle(ipart).type();
    const int nsite1 = particle(ipart).num_sites();
    const int type2 = config.particle(ipart).type();
    const int nsite2 = config.particle(ipart).num_sites();
    ASSERT(type1 == type2, "type1: " << type1 << " != type2 " << type2);
    ASSERT(nsite1 == nsite2, "nsite1: " << nsite1 << " != nsite2 " << nsite2);
  }

  if (add_missing) {
    for (int extra = num_particle_types();
         extra < config.num_particle_types();
         ++extra) {
      add_particle_type(config.type_to_file_name(extra));
    }
    for (int extra = num_particles(); extra < config.num_particles(); ++extra) {
      add_particle_of_type(config.particle(extra).type());
    }
  }
  ASSERT(num_particles() == config.num_particles(),
    "number of particles in config " << num_particles() <<
    " does not match number of given config " << config.num_particles());

  DEBUG("replace positions");
  int part2 = 0;
  for (const int part1 : selection_of_all().particle_indices()) {
    replace_position_(part1, config.particle(part2));
    ++part2;
  }
}

int Configuration::dimension() const { return domain().dimension(); }

void Configuration::synchronize_(const Configuration& config,
    const Select& perturbed) {
  DEBUG(perturbed.str());
  Select sync_sel(perturbed, config.particles_);
  update_positions(sync_sel, true, true);
//  for (int spindex = 0; spindex < perturbed.num_particles(); ++spindex) {
//    const int part_index = perturbed.particle_index(spindex);
//    one_site_select_.set_particle(0, part_index);
//    DEBUG("part_index " << part_index);
//    const Particle& part = config.select_particle(part_index);
//    for (int sindex : perturbed.site_indices(spindex)) {
////      DEBUG("old pos: " << select_particle(part_index).site(sindex).position().str());
////      Site * site = particles_.get_particle(part_index)->get_site(sindex);
////      *site = part.site(sindex);
////      DEBUG("new pos: " << select_particle(part_index).site(sindex).position().str());
////      position_tracker_(part_index, sindex);
////    }
//    *particles_.get_particle(part_index) = part;
//    position_tracker_(part_index);
//  }
}

void Configuration::set_particle_type(const int ptype,
                                      const Select& select) {
  for (int particle_index : select.particle_indices()) {
    Particle * part = particles_.get_particle(particle_index);
    --num_particles_of_type_[part->type()];
    part->set_type(ptype);
    ++num_particles_of_type_[part->type()];
    ASSERT(part->num_sites() == particle_type(ptype).num_sites(),
      "cannot morph into particle with different number of sites");
    for (int isite = 0; isite < part->num_sites(); ++isite) {
      part->get_site(isite)->set_type(particle_type(ptype).site(isite).type());
    }
    for (Select& select : group_selects_) {
      update_selection_(particle_index, &select);
    }
    // HWH doesn't update type-based cell lists, groups, etc.
  }
}

void Configuration::change_volume(const double delta_volume,
    const argtype& args) {
  ASSERT(domain().volume() + delta_volume > 0,
    "delta_volume " << delta_volume << " too large for volume "
    << domain().volume());
  ASSERT(wrap_, "positions must be wrapped before scaling for volume change");
  double factor = 1. + delta_volume/domain_->volume();
  Arguments args_(args);
  const int dimen = args_.key("dimension").dflt("-1").integer();
  if (dimen == -1) {
    factor = std::pow(factor, 1./dimension());
    for (int dim = 0; dim < dimension(); ++dim) {
      domain_->set_side_length(dim, domain_->side_length(dim)*factor);
    }
  } else {
    domain_->set_side_length(dimen, domain_->side_length(dimen)*factor);
  }
  particles_.scale_particle_positions(dimen, factor);
  position_tracker_();
}

}  // namespace feasst
