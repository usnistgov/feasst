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
  ASSERT(partial_configs_.size() == 0, "types cannot be added after groups");
  ASSERT(num_particles() == 0, "types cannot be added after particles");
  particle_types_.add(file_name);
  unique_types_.add(file_name);
}

void Configuration::add_(const Particle particle) {
  Particle part = particle;
  part.update_cell(domain());
  particles_.add(part);
  for (Configuration& partial : partial_configs_) {
    const Group& group = partial.group_;
    if (group.dynamic()) {
      if (group.is_in(part)) {
        std::vector<int> p2f_sitemap, f2p_sitemap;
        Particle filtered = group.remove_sites(part,
                                               &f2p_sitemap,
                                               &p2f_sitemap);
        /// Update cell list
        filtered.update_cell(domain());

        partial.add_(filtered);
        partial.partial_to_full_.push_back(num_particles() - 1);
        partial.full_to_partial_.push_back(partial.num_particles() - 1);
        partial.partial_to_full_site_.push_back(p2f_sitemap);
        partial.full_to_partial_site_.push_back(f2p_sitemap);
      } else {
        partial.full_to_partial_.push_back(-1);
        partial.full_to_partial_site_.push_back(std::vector<int>());
      }
    }
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
  for (Configuration& partial : partial_configs_) {
    const Group& group = partial.group_;
    if (group.dynamic()) {
      std::vector<int>* ftp = &partial.full_to_partial_;
      std::vector<std::vector<int> >* ftps = &partial.full_to_partial_site_;
      const int partial_particle_index = (*ftp)[particle_index];
      if (partial_particle_index >= 0) {
        partial.remove_particle_(partial_particle_index);

        // decrement partial_to_full_
        std::vector<int>* ptf = &partial.partial_to_full_;
        std::vector<std::vector<int> >* ptfs = &partial.partial_to_full_site_;
        for (int index = partial_particle_index;
             index < partial.num_particles();
             ++index) {
          --(*ptf)[index];
        }
        ptf->erase(ptf->begin() + partial_particle_index);
        ptfs->erase(ptfs->begin() + partial_particle_index);

        // decrement full_to_partial_
        for (int index = particle_index; index < num_particles(); ++index) {
          if ((*ftp)[index] != -1) {
            --(*ftp)[index];
          }
        }
      }
      ftp->erase(ftp->begin() + particle_index);
      ftps->erase(ftps->begin() + particle_index);
    }
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

  // do the same for the partial configurations.
  for (Configuration& partial : partial_configs_) {
    const Group& group = partial.group_;
    if (group.dynamic()) {
      const int partial_particle_index =
        partial.full_to_partial_[particle_index];
      if (partial_particle_index >= 0) {
        const Particle part = particles_.particle(particle_index);
        const Particle filtered = group.remove_sites(part);
        partial.replace_position_(partial_particle_index, filtered);
      }
    }
  }
  position_tracker_(particle_index);
}

void Configuration::displace_site_(const int particle_index,
                                   const int site_index,
                                   const Position &displacement) {
  Position pos = displacement;
  pos.add(particle(particle_index).site(site_index).position());
  particles_.replace_position(particle_index, site_index, pos);
  for (Configuration& partial : partial_configs_) {
    const Group& group = partial.group_;
    if (group.dynamic()) {
      const int partial_particle_index =
        partial.full_to_partial_[particle_index];
      if (partial_particle_index >= 0) {
        const int partial_site_index =
          partial.full_to_partial_site_[particle_index][site_index];
        partial.replace_position_(partial_particle_index,
                                  partial_site_index,
                                  pos);
      }
    }
  }
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
  for (Configuration& partial : partial_configs_) {
    const Group& group = partial.group_;
    if (group.dynamic()) {
      const int partial_particle_index =
        partial.full_to_partial_[particle_index];
      if (partial_particle_index >= 0) {
        partial.replace_position_(partial_particle_index,
                                  group.remove_sites(replacement));
      }
    }
  }
  position_tracker_(particle_index);
}

void Configuration::replace_position_(const int particle_index,
                                      const int site_index,
                                      const Position& replacement) {
  particles_.replace_position(particle_index, site_index, replacement);
  for (Configuration& partial : partial_configs_) {
    const Group& group = partial.group_;
    if (group.dynamic()) {
      const int partial_particle_index =
        partial.full_to_partial_[particle_index];
      if (partial_particle_index >= 0) {
        const int partial_site_index =
          partial.full_to_partial_site_[particle_index][site_index];
        partial.replace_position_(partial_particle_index,
                                  partial_site_index,
                                  replacement);
      }
    }
  }
  position_tracker_(particle_index, site_index);
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
  partial_configs_.push_back(*this);
  Configuration * partial = &partial_configs_.back();
  partial->group_ = group;
  partial->particles_.remove(group);
  partial->particle_types_.remove(group);
  partial->unique_types_.remove(group);
  ASSERT(group_.empty(), "main configuration contains no groups");
}

void Configuration::position_tracker_(const int particle_index,
                                      const int site_index) {
  /// update cells
  particles_.update_cell(domain(), particle_index, site_index);
  for (Configuration& partial : partial_configs_) {
    const int partial_particle_index =
      partial.full_to_partial_[particle_index];
    if (partial_particle_index >= 0) {
      // update cell of entire particle if site_index is -1
      if (site_index == -1) {
        partial.particles_.update_cell(partial.domain(),
                                       partial_particle_index);
      } else {
        const int partial_site_index =
          partial.full_to_partial_site_[particle_index][site_index];
        partial.particles_.update_cell(partial.domain(),
                                       partial_particle_index,
                                       partial_site_index);
      }
    }
  }
  /// update neighbors? groups?
}

void Configuration::position_tracker_() {
  for (int index = 0; index < num_particles(); ++index) {
    position_tracker_(index);
  }
}

void Configuration::init_cells(const double min_length,
                               const int partial_index) {
  if (partial_index == -1) {
    domain_.init_cells(min_length);
  } else {
    partial_configs_[partial_index].domain_.init_cells(min_length);
  }
  position_tracker_();
}

/// HWH add check_size .. domain, positions, particles, etc
void Configuration::check_size() const {
  particles_.check_size();
  particle_types_.check_size();
  unique_types_.check_size();
  for (const Configuration& partial : partial_configs_) {
    ASSERT(static_cast<int>(partial.full_to_partial_.size()) == num_particles(),
      "size error");
    partial.check_size();
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
  int random_index = -1;

  // do not select if no particles
  if (num_particles() == 0) {
    selection_.clear();
    return;
  }

  // randomly select a few particles in hopes of finding one of type
  bool term = false;
  int iteration = 0;
  const int max_iteration = 10;
  while (!term && (iteration < max_iteration)) {
    random_index = random_.uniform(0, num_particles() - 1);
    if (particle(random_index).type() == type) {
      term = true;
    }
    ++iteration;
  }

  // if particle of type is not found, enumerate all choices.
  if (iteration == max_iteration) {
    std::vector<int> parts;
    for (int index = 0; index < num_particles(); ++index) {
      if (particle(index).type() == type) {
        parts.push_back(index);
      }
    }

    if (parts.size() == 0) {
      selection_.clear();
      return;
    } else {
      random_index = random_.element(parts);
    }
  }

  select_particle(random_index);
}

void Configuration::set_selection(const Selection selection) {
  std::string id = selection.unique_id();
  ASSERT(id.empty() || id == unique_indices_,
    "If selection is obtained from a configuration, it will have a unique" <<
    "identifier attached to it. The id of the selection and the " <<
    "configuration were fount to be inconsistent, likely due to an " <<
    "expired selection (e.g., particles were removed after selection)");
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

void Configuration::load_coordinates(
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
  for (const std::pair<int, std::vector<int>> pair : selection_.selection()) {
    const int particle_index = pair.first;
    for (const int site_index : pair.second) {
      displace_site_(particle_index, site_index, displacement);
    }
  }
}

}  // namespace feasst
