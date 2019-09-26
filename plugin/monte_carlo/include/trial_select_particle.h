
#ifndef FEASST_MONTE_CARLO_TRIAL_SELECT_PARTICLE_H_
#define FEASST_MONTE_CARLO_TRIAL_SELECT_PARTICLE_H_

#include <vector>
#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_select.h"

namespace feasst {

/// Select a random particle for trial.
class TrialSelectParticle : public TrialSelect {
 public:
  TrialSelectParticle(
    /**
      load_coordinates : load the coordinates into the selection (default: true)

      site : site index to select. If all sites, set to -1 (default).

      ghost : select ghost particles (default: false).
     */
    const argtype& args = argtype()) : TrialSelect(args) {
    class_name_ = "TrialSelectParticle";
    Arguments args_(args);
    args_.dont_check();
    load_coordinates_ = args_.key("load_coordinates").dflt("true").boolean();

    // parse site
    site_ = args_.key("site").dflt("-1").integer();
    if (site_ != -1) {
      mobile_.clear();
      mobile_.add_site(0, site_);
      site_vec_ =  {site_};
    }

    set_ghost(args_.key("ghost").dflt("false").boolean());
  }

  /// Return true if loading coordinates into selection.
  bool load_coordinates() const { return load_coordinates_; }

  /// Add random particle in group index to select.
  /// Return the number of particles to choose from.
  int random_particle(const Configuration& config,
      SelectPosition * select,
      Random * random) {
    ASSERT(group_index() >= 0, "error");
    const int num = config.num_particles(group_index());
    if (num > 0) {
      const int index = random->uniform(0, num - 1);
      const SelectGroup& ran = config.group_select(group_index());
      DEBUG("index " << group_index() << " " << index);
      DEBUG("num " << ran.num_particles());
      bool fast;
      if (site_ == - 1) {
        fast = select->replace_indices(ran.particle_index(index),
                                       ran.site_indices(index));
      } else {
        fast = select->replace_indices(ran.particle_index(index),
                                       site_vec_);
      }
      if (load_coordinates()) {
        if (!fast) select->resize();
        select->load_positions(config.particles());
      }
    } else {
      select->clear();
    }
    return num;
  }

  /// Add random particle in group index to select.
  /// Return the number of particles to choose from.
  void ghost_particle(Configuration * config,
    SelectPosition * select) {
    ASSERT(static_cast<int>(config->ghosts().size()) > particle_type(),
      "type not recognized");
    // if no ghosts, create one
    DEBUG("num ghosts " << config->ghosts()[particle_type()].num_particles());
    if (config->ghosts()[particle_type()].num_particles() == 0) {
      config->add_particle_of_type(particle_type());
      Select add;
      add.add_particle(config->newest_particle(), config->newest_particle_index());
      config->remove_particles(add);
    }
    const Select& ghost = config->ghosts()[particle_type()];
    bool fast;
    // replace indices with the last ghost as may be optimal method available
    // to delete.
    if (site_ == -1) {
      fast = select->replace_indices(ghost.particle_indices().back(),
                                     ghost.site_indices().back());
    } else {
      fast = select->replace_indices(ghost.particle_indices().back(),
                                     site_vec_);
      config->set_selection_physical(ghost, false);
      config->set_selection_physical(*select, true);
    }
    if (load_coordinates()) {
      if (!fast) select->resize();
      select->load_positions(config->particles());
    }
  }

  bool select(const Select& perturbed, System* system, Random * random) override {
    if (is_ghost()) {
      ghost_particle(system->get_configuration(), &mobile_);
      set_probability(1.);
    } else {
      const int num = random_particle(system->configuration(), &mobile_, random);
      if (num <= 0) return false;
      set_probability(1./static_cast<double>(num));
    }
    mobile_.remove_unphysical_sites(system->configuration());
    mobile_original_ = mobile_;
    return true;
  }

  std::shared_ptr<TrialSelect> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialSelectParticle(std::istream& istr);
  virtual ~TrialSelectParticle() {}

 protected:
  void serialize_trial_select_particle_(std::ostream& ostr) const;

 private:
  bool load_coordinates_;
  int site_;
  std::vector<int> site_vec_;
};

inline std::shared_ptr<TrialSelectParticle> MakeTrialSelectParticle(
    const argtype &args = argtype()) {
  return std::make_shared<TrialSelectParticle>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_SELECT_PARTICLE_H_
