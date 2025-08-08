
#ifndef FEASST_MONTE_CARLO_TRIAL_SELECT_PARTICLE_H_
#define FEASST_MONTE_CARLO_TRIAL_SELECT_PARTICLE_H_

#include <map>
#include <string>
#include <vector>
#include <memory>
#include "monte_carlo/include/trial_select.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/// Select a random particle for Trial.
class TrialSelectParticle : public TrialSelect {
 public:
  //@{
  /*
    Users typically do not set the following values:
    - load_coordinates: load the coordinates into the selection (default: true)
    - ghost: select ghost particles (default: false).
    - half_ghost: select ghost particles half of the time (default: false).
    - exclude_perturbed: if true, exclude perturbed particle (default: false)
   */
  /** @name Arguments
    - site: site name to select. If all sites, set to -1 (default).
    - min_particles: do not select if less than min number of particles.
      If -1, no constraint (default: -1).
      Note that this is the number of particles before the Trial.
    - max_particles: do not select if more than max number of particles.
      If -1, no constraint (default: -1).
      Note that this is the number of particles before the Trial.
    - TrialSelect arguments.
   */
  explicit TrialSelectParticle(argtype args = argtype());
  explicit TrialSelectParticle(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void precompute(System * system) override;

  /// Return true if loading coordinates into selection.
  bool load_coordinates() const { return load_coordinates_; }

  /// Return site.
  int site() const { return site_; }

  /// Add random particle in group index to select.
  /// Return the number of particles to choose from.
  int random_particle(const Configuration& config,
    /// Exclude from selection.
    const Select * exclude,
    Select * select,
    Random * random);

  /// Same as above, but no exclusions.
  int random_particle(const Configuration& config,
    Select * select,
    Random * random) { return random_particle(config, NULL, select, random); }

  /// Select a ghost particle.
  void ghost_particle(Configuration * config,
    /// Exclude from selection.
    const Select* exclude,
    Select * select);

  /// Same as above, but no exclusions.
  void ghost_particle(Configuration * config, Select * select) {
    ghost_particle(config, NULL, select); }

  bool select(const Select& perturbed,
              System* system,
              Random * random,
              TrialSelect * previous_select) override;

  /// Select a particular particle by index.
  /// Note that this index ignores ghost particles.
  void select_particle(const int particle_index, const Configuration& config);

  std::shared_ptr<TrialSelect> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialSelectParticle(std::istream& istr);
  virtual ~TrialSelectParticle() {}

  //@}
 protected:
  void serialize_trial_select_particle_(std::ostream& ostr) const;

 private:
  bool load_coordinates_;
  bool exclude_perturbed_;
  int site_ = -123; // avoid uninitialized warning but requires precompute
  std::string site_name_;
  int min_particles_;
  int max_particles_;
  std::vector<int> site_vec_;
  bool half_ghost_;

  // compute the number of excluded particles which are candidates for selection
  int num_excluded_(const Configuration& config, const Select * exclude);
};

inline std::shared_ptr<TrialSelectParticle> MakeTrialSelectParticle(
    argtype args = argtype()) {
  return std::make_shared<TrialSelectParticle>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_SELECT_PARTICLE_H_
