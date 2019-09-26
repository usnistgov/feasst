
#ifndef FEASST_SYSTEM_SELECT_LIST_H_
#define FEASST_SYSTEM_SELECT_LIST_H_

#include "configuration/include/configuration.h"
#include "math/include/utils_math.h" // move to cpp

namespace feasst {

/**
  Similar to Select, except includes configurations as argument.
    Selects are for modifications of a subset of the configuration.
    Note that the requested selection may not be possible.
    In this case the selection remains empty.
 */
// HWH: contain SelectPosition instead of inheriting from it to hide selection interface
// HWH: consider depreciating/moving SelectList
class SelectList : public SelectPosition {
 public:
  SelectList() {}

  SelectList& particle(
      // HWH document this index.
      const int index,
      const Configuration& config,
      // By default, if group_index is 0, consider all particles
      const int group_index = 0) {
    // HWH optimize this
    clear();
    const SelectGroup& select = config.group_selects()[group_index];
    ASSERT(index < select.num_particles(), "error");
    add_particle(select.particle_index(index), select.site_indices(index));
    resize();
    load_positions(config.particles());
    return *this;
  }

  // fast replace of single particle
  void replace_particle(const Select& replacement,
      const int sp_index, // selection particle index
      const Configuration& config) {
    bool fast = replace_indices(replacement.particle_index(sp_index),
                                replacement.site_indices(sp_index));
    if (!fast) resize();
    load_positions(config.particles());
  }

  /// Select the last particle that was added to configuration.
  // HWH depreciate
  void last_particle_added(const Configuration * config) {
    // HWH optimize for add delete trial
    clear();
    add_particle(config->newest_particle(), config->newest_particle_index());
    resize();
    load_positions(config->particles());
  }

  /// Select sites.
  void select_sites(const Configuration& config,
                    const int particle_index,
                    const std::vector<int> site_indices) {
    clear();
    // filter the particle index.
    const int index = config.selection_of_all().particle_index(particle_index);
    add_sites(index, site_indices);
    resize();
    load_positions(config.particles());
  }

  /// Select a site.
  void select_site(const Configuration& config,
                   const int particle_index,
                   const int site_index) {
    select_sites(config, particle_index, {site_index});
  }

  /// Note: this method is relatively unoptimized compared with other options.
  virtual void add(const Configuration& config, const Group& group) {
    for (int index = 0; index < config.num_particles(); ++index) {
      const Particle& part = config.particle(index);
      if (group.is_in(part)) {
        Particle filtered = part;
        group.remove_sites(&filtered);
        add_particle(filtered, index);
      }
    }
  }

  /// Return the selected particle.
  const Particle& particle(const Configuration& config) {
    ASSERT(num_particles() == 1, "requires one selected particle");
    return config.select_particle(particle_index(0));
  }

  /// Return the bond between two sites in particle use for model parameters
  const Bond& bond(
    /// The particle index is the one in the selection, not the configuration
    const int particle,
    const int site1, const int site2, const Configuration& config) const {
    const int particle_index = particle_indices()[particle];
    const Particle& part = config.select_particle(particle_index);
    const int particle_type = part.type();
    const Particle& type = config.particle_types().particle(particle_type);
    const int bond_type = type.bond(site1, site2).type();
    const Particle& unique = config.unique_types().particle(particle_type);
    return unique.bond(bond_type);
  }

  /// Select the current positions of given selection
  void store(const Select& select, const Configuration& config) {
    clear();
    Select::add(select);
    resize();
    load_positions(config.particles());
  }

  /// Remove unphysical sites from selection
  void remove_unphysical_sites(const Configuration& config) {
    Select unphysical;
    for (int sp_index = 0;
         sp_index < static_cast<int>(particle_indices().size());
         ++sp_index) {
      const int p_index = particle_indices()[sp_index];
    //for (const int p_index : particle_indices()) {
      std::vector<int> sites;
      for (const int s_index : site_indices(sp_index)) {
        if (!config.select_particle(p_index).site(s_index).is_physical()) {
          sites.push_back(s_index);
        }
      }
      if (sites.size() > 0) {
        unphysical.add_sites(p_index, sites);
      }
    }
    if (unphysical.num_particles() > 0) {
      remove(unphysical);
      resize();
      load_positions(config.particles());
    }
  }

  void serialize(std::ostream& ostr) const override {
    SelectPosition::serialize(ostr); }
  SelectList(std::istream& istr) : SelectPosition(istr) {}
  virtual ~SelectList() {}
};

}  // namespace feasst

#endif  // FEASST_SYSTEM_SELECT_LIST_H_
