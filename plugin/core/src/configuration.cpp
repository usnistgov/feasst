#include <fstream>
#include <sstream>
#include "core/include/configuration.h"
#include "core/include/file_xyz.h"
#include "core/include/utils.h"
#include "core/include/utils_math.h"
#include "core/include/debug.h"
#include "core/include/constants.h"

namespace feasst {

Configuration::Configuration(const argtype& args) {
  particle_types_.unique_particles();
  unique_types_.unique_types();
  reset_unique_indices_();
  add(Group()); // add empty group which represents all particles
  args_.init(args);

  // parse domain
  if (args_.key("cubic_box_length").used()) {
    set_domain(Domain().set_cubic(args_.dble()));
  }

  // parse types
  std::string start("particle_type");
  // if only one particle type, allow drop the subscript
  if (args_.key(start).used()) {
    add_particle_type(args_.str());
  } else {
    std::stringstream key;
    int type = num_particle_types();
    key << start << type;
    while (args_.key(key.str()).used()) {
      add_particle_type(args_.str());
      ++type;
      ASSERT(type < 1e8, "type(" << type << ") is very high. Infinite loop?");
      key.str("");
      key << start << type;
    }
  }

  if (args_.key("init_cells").used()) {
    const double min_length = args_.dble();
    int group_index = args_.key("cell_group").dflt("0").integer();
    init_cells(min_length, group_index);
  }
}

void Configuration::reset_unique_indices_() {
  unique_indices_ = random_.alpha_numeric();
}

void Configuration::add_particle_type(const std::string file_name) {
  DEBUG("adding type");
  ASSERT(num_particles() == 0, "types cannot be added after particles");
  particle_types_.add(file_name);
  unique_types_.add(file_name);
  ghosts_.push_back(SelectGroup());
  ASSERT(ghosts_.back().group().is_empty(), "");
  ASSERT(!find_in_list(file_name, type_to_file_),
    "file_name(" << file_name << ") already provided.");
  type_to_file_.push_back(file_name);
  num_particles_of_type_.push_back(0);
}

void Configuration::add_(const Particle particle) {
  Particle part = particle;
  particles_.add(part);
  for (SelectGroup& select : group_selects_) {
    add_to_selection_(num_particles() - 1, &select);
  }
  position_tracker_(num_particles() - 1);
}

void Configuration::add_particle_of_type(const int type) {
  ASSERT(type < num_particle_types(), "type(" << type << ") is not allowed "
    << "when there are only " << num_particle_types() << " particle types");
  // decide whether to add a ghost particle or a new one
  if (ghosts_[type].num_particles() == 0) {
    Particle part = particle_types_.particle(type);
    part.erase_bonds();
    add_(part);
    newest_particle_index_ = num_particles() - 1;
  } else {
    const int index = ghosts_[type].particle_index(0);
    ghosts_[type].remove_particle(index);
    for (SelectGroup& select : group_selects_) {
      add_to_selection_(index, &select);
    }
    newest_particle_index_ = index;
  }
  ++num_particles_of_type_[type];
}

void Configuration::remove_particle_(const int particle_index) {
  reset_unique_indices_();
  const int type = particles_.particle(particle_index).type();
  --num_particles_of_type_[type];
  ghosts_[type].add_particle(select_particle(particle_index), particle_index);
  DEBUG("type " << type);
  DEBUG("particle index " << particle_index);
  DEBUG("num particles " << num_particles());
  for (SelectGroup& select : group_selects_) {
    select.remove_particle(particle_index);
  }

  // remove particle_index from cell
  // note: somewhat derivative of position_tracker
  for (int cell_index = 0;
       cell_index < static_cast<int>(domain().cells().size());
       ++cell_index) {
    const Cells& cells = domain().cells()[cell_index];
    const int group_index = feasst::round(cells.property("group"));
    const Particle& part = select_particle(particle_index);
    const Group& group = group_selects_[group_index].group();
    if (group.is_in(part)) {
      for (int site_index = 0; site_index < part.num_sites(); ++site_index) {
        const Site& site = part.site(site_index);
        if (group.is_in(site)) {
          const std::string name = cells.label();
          double value;
          if (site.properties().value(name, &value)) {
            const int cell_old = feasst::round(value);
            Select select;
            select.add_site(particle_index, site_index);
            domain_.remove_from_cell_list(cell_index, select, cell_old);
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
  // first, obtain the particle position
  Position position = select_particle(particle_index).position();
  position.add(displacement);

  // apply boundary conditions on the new position
  const Position pbc_shift = domain().shift(position);

  // obtain the potentially new displacement in light of the boundaries
  Position new_displacement(displacement);
  new_displacement.add(pbc_shift);

  // add the new displacement to the particle position
  particles_.displace(particle_index, new_displacement);

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

void Configuration::add(Group group, std::string name) {
  ASSERT(group.is_empty() || particle_types_.num() != 0,
    "add groups after particle types");
  if (name == "") {
    std::stringstream ss;
    ss << num_groups();
    name.assign(ss.str());
  }
  SelectGroup group_select;
  group.add_property(name, 0.);
  group_select.set_group(group);
  init_selection_(&group_select);
  group_selects_.push_back(group_select);
}

void Configuration::position_tracker_(const int particle_index,
                                      const int site_index) {
  ASSERT(site_index >= 0, "index error");
  DEBUG("update cells");
  for (int cell_index = 0;
       cell_index < static_cast<int>(domain().cells().size());
       ++cell_index) {
    DEBUG("cell index " << cell_index);
    const Cells& cells = domain().cells()[cell_index];
    DEBUG("group " << cells.properties().value("group"));
    const int group_index = feasst::round(cells.property("group"));
    ASSERT(group_index >= 0, "error");
    const Particle& part = select_particle(particle_index);
    const Group& group = group_selects_[group_index].group();
    if (group.is_in(part)) {
      const Site& site = part.site(site_index);
      if (group.is_in(site)) {
        const int cell_new = domain().cell_id(site.position(), cells);
        Select select;
        select.add_site(particle_index, site_index);
        const std::string name = cells.label();
        double value;
        if (site.properties().value(name, &value)) {
          const int cell_old = feasst::round(value);
          DEBUG("index " << particle_index << " " << site_index);
          DEBUG("new cell " << cell_new << " old cell " << cell_old);
          DEBUG("before new cell set: " << particles_.particle(particle_index).site(site_index).property("cell0"));
          particles_.set_site_property(name, cell_new, particle_index, site_index);
          domain_.update_cell_list(cell_index, select, cell_new, cell_old);
        } else {
          particles_.add_site_property(name, cell_new, particle_index, site_index);
          DEBUG("adding to cell list " << cell_index << " cllnw " << cell_new << " si " << site_index);
          domain_.add_to_cell_list(cell_index, select, cell_new);
        }
      }
    }
  }

  DEBUG("update selection");
  for (const SelectGroup& select : group_selects_) {
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
  for (int index = 0; index < num_particles(); ++index) {
    position_tracker_(index);
  }
}

void Configuration::init_cells(const double min_length,
                               const int group_index) {
  domain_.init_cells(min_length, group_index);
  position_tracker_();
}

/// HWH add check .. domain, positions, particles, etc
void Configuration::check() const {
  particles_.check();
  particle_types_.check();
  unique_types_.check();

  ASSERT(particle_types_.num() == num_particle_types(), "er");
  ASSERT(unique_types_.num_sites() == num_site_types(), "er");

  // check that the first group is all particles in the configuration.
  ASSERT(static_cast<int>(group_selects_[0].num_particles()) == num_particles(),
    "size error");
  for (int index = 0; index < num_particles(); ++index) {
    ASSERT(static_cast<int>(group_selects_[0].site_indices(index).size()) ==
      particle(index).num_sites(), "size error");
  }

  // check that number of particles in cell list is number in selection.
  for (const Cells& cells : domain().cells()) {
    const int group_index = feasst::round(cells.property("group"));
    const Select& select = group_selects_[group_index];
    ASSERT(select.num_sites() == cells.num_sites(),
      "sites in group(" << select.num_sites() << ") is not equal to sites " <<
      "in cell(" << cells.num_sites() << ")");
  }

  // check number of particle types
  ASSERT(
    static_cast<int>(num_particles_of_type_.size()) == num_particle_types(),
    "size error"
  );

  model_params().check();
}

void Configuration::check_id_(const Select* select) const {
  check_id_(select->unique_id());
}
void Configuration::check_id_(const Select& select) const {
  check_id_(select.unique_id());
}
void Configuration::check_id_(const std::string id) const {
  ASSERT(id.empty() || id == unique_indices_,
    "If selection is obtained from a configuration, it will have a unique" <<
    "identifier attached to it. The id of the selection and the " <<
    "configuration were fount to be inconsistent, likely due to an " <<
    "expired selection (e.g., particles were removed after selection)");
}

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
                                      SelectGroup * select) const {
  const Particle& part = select_particle(particle_index);
  const Group group = select->group();
  if (group.is_in(part)) {
    select->add_particle(particle_index, group.site_indices(part));
  }
}

void Configuration::init_selection_(SelectGroup * group_select) const {
  for (int part_index = 0; part_index < num_particles(); ++part_index) {
    add_to_selection_(part_index, group_select);
  }
}

void Configuration::update_positions(const SelectPosition& select) {
  check_id_(select);
  int pindex = 0;
  for (int particle_index : select.particle_indices()) {
    replace_position_(particle_index, select.particle_positions()[pindex]);
    int sindex = 0;
    for (int site_index : select.site_indices(pindex)) {
      replace_properties_(particle_index,
                          site_index,
                          select.site_properties()[pindex][sindex],
                          {"cell"});
      replace_position_(particle_index,
                        site_index,
                        select.site_positions()[pindex][sindex]);
      ++sindex;
    }
    ++pindex;
  }
}

int Configuration::particle_type_to_group(const int particle_type) {
  int index;
  if (num_particle_types() == 1) {
    return 0;
  }
  if (!find_in_list(particle_type, group_store_particle_type_, &index)) {
    index = num_groups();
    add(Group().add_particle_type(particle_type));
    group_store_particle_type_.push_back(particle_type);
    group_store_group_index_.push_back(index);
  }
  return index;
}

int Configuration::num_particles(const int group) const {
  TRACE("pn " << particles_.num());
  TRACE("gn " << num_ghosts_());
  TRACE("group " << group);
  if (group == 0) {
    TRACE("here " << particles_.num() - num_ghosts_());
    return particles_.num() - num_ghosts_();
  }
  TRACE("here " << group_selects_[group].num_particles());
  return group_selects_[group].num_particles();
}

int Configuration::num_ghosts_() const {
  int num = 0;
  for (const SelectGroup& select : ghosts_) {
    num += select.num_particles();
  }
  return num;
}

void Configuration::revive(const SelectPosition& selection) {
  for (int particle_index : selection.particle_indices()) {
    const Particle& part = select_particle(particle_index);
    const int type = part.type();
    ghosts_[type].remove_last_particle();
    for (SelectGroup& select : group_selects_) {
      add_to_selection_(particle_index, &select);
    }
  }
}

const Particle Configuration::particle(const int index,
                                       const int group) const {
  const SelectGroup& select_group = group_selects_[group];
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

}  // namespace feasst
