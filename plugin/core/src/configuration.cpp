#include <fstream>
#include <sstream>
#include <utility>  // pair
#include "core/include/configuration.h"
#include "core/include/file_xyz.h"
#include "core/include/utils.h"
#include "core/include/debug.h"
#include "core/include/constants.h"

namespace feasst {

Configuration::Configuration() {
  particle_types_.unique_particles();
  unique_types_.unique_types();
  reset_unique_indices_();
}

void Configuration::reset_unique_indices_() {
  unique_indices_ = random_.alpha_numeric();
  selection_.set_unique_id(unique_indices_);
}

void Configuration::add_particle_type(const char* file_name) {
  DEBUG("adding type");
  ASSERT(num_particles() == 0, "types cannot be added after particles");
  particle_types_.add(file_name);
  unique_types_.add(file_name);
}

void Configuration::add_(const Particle particle) {
  Particle part = particle;
  part.update_cell(domain());
  particles_.add(part);
  select_all_.add_last_particle(particles_);
  for (GroupSelection& select : group_selects_) {
    add_to_selection_(num_particles() - 1, &select);
  }
}

void Configuration::add_particle(const int type) {
  ASSERT(type < num_particle_types(), "type(" << type << ") is not allowed "
    << "when there are only " << num_particle_types() << " particle types");
  Particle part = particle_types_.particle(type);
  part.erase_bonds();
  add_(part);
}

void Configuration::remove_particle_(const int particle_index) {
  reset_unique_indices_();
  particles_.remove(particle_index);
  select_all_.remove_particle(particle_index);
  DEBUG("particle index " << particle_index);
  DEBUG("num particles " << num_particles());
  for (GroupSelection& select : group_selects_) {
    select.remove_particle(particle_index);
  }
}

void Configuration::remove_selected_particles() {
  ASSERT(selection_.num_particles() > 0, "no selection");
  // loop through selection backwards.
  // HWH: sort selection
  ASSERT(selection_.num_particles() == 1, "implement sort");
  for (int index = selection_.num_particles()  - 1;
       index >= 0;
       --index) {
    remove_particle_(selection_.particle_index(index));
  }
  selection_.clear();
}

void Configuration::remove_selected_particle() {
  ASSERT(selection_.num_particles() == 1, "assumes 1 particle selected");
  remove_selected_particles();
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

void Configuration::displace_selected_particles(const Position &displacement) {
  ASSERT(selection_.num_particles() > 0, "no selection");
  for (int particle_index : selection_.particle_indices()) {
    displace_particle_(particle_index, displacement);
  }
}

void Configuration::displace_selected_particle(const Position &displacement) {
  ASSERT(selection_.num_particles() == 1, "assumes 1 particle selected");
  displace_selected_particles(displacement);
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

void Configuration::replace_selected_particle_position(
    const Particle& replacement) {
  ASSERT(selection_.num_particles() == 1, "requires one particle selected");
  replace_position_(selection_.particle_index(0), replacement);
}

void Configuration::replace_position_of_last(const Particle& replacement) {
  replace_position_(num_particles() - 1, replacement);
}

void Configuration::default_configuration() {
  particle_types_ = Particles().unique_particles().default_particles();
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

void Configuration::add(const Group group) {
  ASSERT(particle_types_.num() != 0, "add groups after particle types");
  GroupSelection group_select;
  group_select.set_group(group);
  init_selection_(&group_select);
  group_selects_.push_back(group_select);
}

void Configuration::position_tracker_(const int particle_index,
                                      const int site_index) {
  /// update cells
  particles_.update_cell(domain(), particle_index, site_index);
  /// HWH update selection?

  /// update neighbors? groups?
}

void Configuration::position_tracker_() {
  for (int index = 0; index < num_particles(); ++index) {
    position_tracker_(index);
  }
}

void Configuration::init_cells(const double min_length) {
  domain_.init_cells(min_length);
  position_tracker_();
}

/// HWH add check_size .. domain, positions, particles, etc
void Configuration::check_size() const {
  particles_.check_size();
  particle_types_.check_size();
  unique_types_.check_size();

  // check that select_all_ is all particles in the configuration.
  ASSERT(static_cast<int>(select_all_.num_particles()) == particles_.num(),
    "size error");
  for (int index = 0; index < particles_.num(); ++index) {
    ASSERT(static_cast<int>(select_all_.site_indices(index).size()) ==
      particle(index).num_sites(), "size error");
  }
}

const Particle& Configuration::selected_particle() const {
  ASSERT(selection_.num_particles() == 1, "requires one particle selected");
  return particle(selection_.particle_index(0));
}

void Configuration::select_particle(const int particle_index) {
  selection_.clear();
  ASSERT(particle_index < num_particles(), "size error");
  selection_.add_particle(particles_, particle_index);
}

void Configuration::select_site(const int particle_index, const int site_index) {
  selection_.clear();
  selection_.add_site(particle_index, site_index);
}

void Configuration::select_last_particle() {
  selection_.clear();
  select_particle(num_particles() - 1);
}

void Configuration::select_random_particle_of_type(const int type) {
  ASSERT(type < num_particle_types(),
    "particle type(" << type << ") doesn't exist");

  // do not select if no particles
  if (num_particles() == 0) {
    selection_.clear();
    return;
  }

  const int group_index = particle_type_to_group_(type);
  select_random_particle_of_group(group_index);
}

void Configuration::check_id_(const Selection& select) const {
  std::string id = select.unique_id();
  ASSERT(id.empty() || id == unique_indices_,
    "If selection is obtained from a configuration, it will have a unique" <<
    "identifier attached to it. The id of the selection and the " <<
    "configuration were fount to be inconsistent, likely due to an " <<
    "expired selection (e.g., particles were removed after selection)");
}

void Configuration::set_selection(const Selection selection) {
  check_id_(selection);
  selection_ = selection;
}

void Configuration::select_random_particle() {
  selection_.clear();
  if (num_particles() > 0) {
    selection_.add_random_particle(particles_);
  }
}

void Configuration::select_all() {
  selection_.clear();
  for (int index = 0; index < num_particles(); ++index) {
    selection_.add_particle(particles_, index);
  }
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

void Configuration::select(const Group& group) {
  for (int index = 0; index < num_particles(); ++index) {
    const Particle& part = particle(index);
    if (group.is_in(part)) {
      Particle filtered = group.remove_sites(part);
      selection_.add_particle(filtered, index);
    }
  }
}

void Configuration::select_random_particle(const Group& group) {
  select(group);
  selection_.reduce_to_random_particle();
}

void Configuration::displace_selection(const Position &displacement) {
  for (const int& particle_index : selection_.particle_indices()) {
    for (const int& site_index : selection_.site_indices(particle_index)) {
      displace_site_(particle_index, site_index, displacement);
    }
  }
}

void Configuration::add_to_selection_(const int particle_index,
                                      GroupSelection * select) const {
  const Particle& part = particle(particle_index);
  const Group group = select->group();
  if (group.is_in(part)) {
    select->add_particle(particle_index, group.site_indices(part));
  }
}

void Configuration::init_selection_(GroupSelection * group_select) const {
  for (int part_index = 0; part_index < num_particles(); ++part_index) {
    add_to_selection_(part_index, group_select);
  }
}

void Configuration::select_random_particle_of_group(const int group_index) {
  ASSERT(group_index < num_groups(), "size error");
  if (group_index == -1) {
    select_random_particle();
  } else {
    selection_ = group_selects_[group_index].random_particle();
  }
  const int num_particles = selection_.num_particles();
  ASSERT(num_particles == 0 || num_particles == 1,
    "only one particle should be selected");
}

void Configuration::update_positions(const PositionSelection& select) {
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

int Configuration::particle_type_to_group_(const int particle_type) {
  int index;
  if (num_particle_types() == 1) {
    return -1;
  }
  if (!find_in_list(particle_type, group_store_particle_type_, &index)) {
    index = num_groups();
    add(Group().add_particle_type(particle_type));
    group_store_particle_type_.push_back(particle_type);
    group_store_group_index_.push_back(index);
  }
  return index;
}

}  // namespace feasst
