#include <fstream>
#include <sstream>
#include <utility>  // pair
#include "core/include/configuration.h"
#include "core/include/file_xyz.h"
#include "core/include/utils.h"
#include "core/include/utils_math.h"
#include "core/include/debug.h"
#include "core/include/constants.h"

namespace feasst {

Configuration::Configuration() {
  particle_types_.unique_particles();
  unique_types_.unique_types();
  reset_unique_indices_();
  add(Group()); // add empty group which represents all particles
}

void Configuration::reset_unique_indices_() {
  unique_indices_ = random_.alpha_numeric();
}

void Configuration::add_particle_type(const char* file_name) {
  DEBUG("adding type");
  ASSERT(num_particles() == 0, "types cannot be added after particles");
  particle_types_.add(file_name);
  unique_types_.add(file_name);
  ghosts_.push_back(SelectGroup());
  ASSERT(ghosts_.back().group().is_empty(), "");
}

void Configuration::add_(const Particle particle) {
  Particle part = particle;
  particles_.add(part);
  for (SelectGroup& select : group_selects_) {
    add_to_selection_(num_particles() - 1, &select);
  }
  position_tracker_(num_particles() - 1);
}

void Configuration::add_particle(const int type) {
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
}

void Configuration::remove_particle_(const int particle_index) {
  reset_unique_indices_();
  const int type = particles_.particle(particle_index).type();
  ghosts_[type].add_particle(particles_, particle_index);
  DEBUG("type " << type);
  DEBUG("particle index " << particle_index);
  DEBUG("num particles " << num_particles());
  for (SelectGroup& select : group_selects_) {
    select.remove_particle(particle_index);
  }
}

void Configuration::remove_particles(const Select& selection) {
  ASSERT(selection.num_particles() > 0, "no selection");
  // loop through selection backwards.
  // HWH: sort selection
  ASSERT(selection.num_particles() == 1, "implement sort");
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
  Position position = particle(particle_index).position();
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
  pos.add(particle(particle_index).site(site_index).position());
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
  ASSERT(particles_.particle(particle_index).num_sites() ==
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

void Configuration::default_configuration() {
  domain_.set_cubic(5.);
  add_particle_type("../forcefield/data.atom");
  add_particle(0);
  add_particle(0);
  Position position = particle(1).site(0).position();
  position.set_coord(0, 1.25);
  Particle part1 = particle(1);
  part1.displace(position);
  replace_position_(1, part1);
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
  /// update cells
  for (const Cells& cells : domain().cells()) {
    DEBUG("group " << cells.properties().value("group"));
    const int group_index = feasst::round(cells.property("group"));
    ASSERT(group_index >= 0, "error");
    if (group_index == 0) {
      // update all particles and sites
      particles_.update_cell(cells, domain(), particle_index, site_index);
    } else {
      const Group& group = group_selects_[group_index].group();
      const Particle& part = particle(particle_index);
      if (group.is_in(part)) {
        for (int gsite : group.site_indices(part)) {
          particles_.update_cell(cells, domain(), particle_index, gsite);
        }
      }
    }
  }
  /// HWH update selection? neighbors?
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

/// HWH add check_size .. domain, positions, particles, etc
void Configuration::check_size() const {
  particles_.check_size();
  particle_types_.check_size();
  unique_types_.check_size();

  // check that the first group is all particles in the configuration.
  ASSERT(static_cast<int>(group_selects_[0].num_particles()) == num_particles(),
    "size error");
  for (int index = 0; index < num_particles(); ++index) {
    ASSERT(static_cast<int>(group_selects_[0].site_indices(index).size()) ==
      particle(index).num_sites(), "size error");
  }
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
  for (int part_index = 0;
       part_index < num_particles();
       ++part_index) {
    int num_site = 0;
    Particle part = particle(part_index);
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
  const Particle& part = particles_.particle(particle_index);
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

}  // namespace feasst
