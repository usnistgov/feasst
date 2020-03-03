
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
  /**
    args:
    - load_coordinates: load the coordinates into the selection (default: true)
    - site: site index to select. If all sites, set to -1 (default).
    - ghost: select ghost particles (default: false).
   */
  explicit TrialSelectParticle(const argtype& args = argtype());

  /// Return true if loading coordinates into selection.
  bool load_coordinates() const { return load_coordinates_; }

  /// Add random particle in group index to select.
  /// Return the number of particles to choose from.
  int random_particle(const Configuration& config,
      Select * select,
      Random * random);

  /// Select a ghost particle.
  void ghost_particle(Configuration * config,
    Select * select);

  bool select(const Select& perturbed,
              System* system,
              Random * random) override;

  /// Select a particular particle by index.
  /// Note that this index ignores ghost particles.
  void select_particle(const int particle_index, const Configuration& config) {
    mobile_.particle(particle_index, config, 0); }

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
