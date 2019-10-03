
#ifndef FEASST_MONTE_CARLO_TRIAL_SELECT_BOND_H_
#define FEASST_MONTE_CARLO_TRIAL_SELECT_BOND_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_select.h"

namespace feasst {

/**
  A random particle of given type is selected if previously perturbed sites
    are not available.
  Select a single bond from given anchor to mobile sites.
 */
class TrialSelectBond : public TrialSelect {
 public:
  TrialSelectBond(
    /**
      mobile_site : index of the mobile site.
      anchor_site : index of the anchor site.
     */
    const argtype& args = argtype());

  int anchor_site() const { return anchor_site_; }
  int mobile_site() const { return mobile_site_; }

  // bond_length is added as a property
  // mobile and anchor are sized
  void precompute(System * system) override {
    TrialSelect::precompute(system);
    const Particle& part = system->configuration().particle_types().particle(particle_type());
    const int bond_type = part.bond(mobile_site_, anchor_site_).type();
    const Bond& bond = system->configuration().unique_types().particle(particle_type()).bond(bond_type);
    add_property("bond_length", bond.property("length"));
    anchor_.clear();
    anchor_.add_site(0, anchor_site_);
    mobile_.clear();
    mobile_.add_site(0, mobile_site_);
  }

  bool select(const Select& perturbed, System * system, Random * random) override {
    Configuration * config = system->get_configuration();
    int particle_index = -1;
    if (perturbed.num_sites() > 0) {
      particle_index = perturbed.particle_indices().back();
      set_probability(1.);
    } else {
      // select random particle of correct type
      const int group_index = config->particle_type_to_group(particle_type());
      const int num = config->num_particles(group_index);
      if (num <= 0) return false;
      const int index = random->uniform(0, num - 1);
      const SelectGroup& select = config->group_select(group_index);
      particle_index = select.particle_index(index);
      set_probability(1./static_cast<double>(num));
    }
    mobile_.set_particle(0, particle_index);
    anchor_.set_particle(0, particle_index);
    mobile_.load_positions(config->particles());
    DEBUG("mobile: " << mobile_.str());
    DEBUG("anchor: " << anchor_.str());
    mobile_original_ = mobile_;
    return true;
  }

  std::shared_ptr<TrialSelect> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialSelectBond(std::istream& istr);
  virtual ~TrialSelectBond() {}

 protected:
  void serialize_trial_select_bond_(std::ostream& ostr) const;

 private:
  int mobile_site_;
  int anchor_site_;
};

inline std::shared_ptr<TrialSelectBond> MakeTrialSelectBond(
    const argtype &args = argtype()) {
  return std::make_shared<TrialSelectBond>(args);
}

}  // namespace feasst

#endif  // FEASST_MONTE_CARLO_TRIAL_SELECT_BOND_H_
