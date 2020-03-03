
#ifndef FEASST_CHAIN_TRIAL_SELECT_SITE_OF_TYPE_H_
#define FEASST_CHAIN_TRIAL_SELECT_SITE_OF_TYPE_H_

#include "monte_carlo/include/trial_select.h"

namespace feasst {

/// Select a random site of a given type.
class TrialSelectSiteOfType : public TrialSelect {
 public:
  /**
    args:
    - site_type: type of site to select.
   */
  TrialSelectSiteOfType(
      const argtype& args = argtype()) : TrialSelect(args) {
    class_name_ = "TrialSelectSiteOfType";
    Arguments args_(args);
    args_.dont_check();
    site_type_ = args_.key("site_type").integer();
  }

  int site_type() const { return site_type_; }

  /// Select a random site of given type in randomly selected particle.
  int random_site_in_particle(
      const Configuration& config,
      Select * select,
      Random * random) {
    ASSERT(config.num_site_types() > site_type(),
      "site_type: " << site_type() << " is not present in system.");
    DEBUG("group_index " << group_index());
    const Select& group = config.group_select(group_index());
    const int num_particles = group.num_particles();
    if (num_particles == 0) {
      DEBUG("no particles");
      return 0;
    }
    const int pindex = random->uniform(0, num_particles - 1);
    const int particle_index = group.particle_index(pindex);
    DEBUG("particle_index " << particle_index);
    const Particle& part = config.select_particle(particle_index);
    const int num_sites_of_type = part.num_sites_of_type(site_type());
    if (num_sites_of_type == 0) {
      DEBUG("no sites of type");
      return 0;
    }
    const int sindex = random->uniform(0, num_sites_of_type - 1);

    // loop through each site in particle to select the one of type
    int count = -1;
    int site_index = -1;
    bool terminate = false;
    while (!terminate) {
      ++site_index;
      if (part.site(site_index).type() == site_type()) {
        ++count;
        if (count == sindex) {
          terminate = true;
        }
      }
      ASSERT(site_index < part.num_sites(), "infinite loop");
    }
    if (select->num_particles() != 1 || select->num_sites(0) != 1) {
      select->clear();
      select->add_site(0, 0);
    }
    select->set_particle(0, particle_index);
    select->set_site(0, 0, site_index);
    DEBUG("selected: " << select->str());
    return num_sites_of_type;
  }

  bool select(const Select& perturbed, System* system, Random * random) override {
    const int num = random_site_in_particle(
      system->configuration(),
      &mobile_,
      random);
    if (num <= 0) return false;
    set_probability(1./static_cast<double>(num));
    mobile_original_ = mobile_;
    return true;
  }

  std::shared_ptr<TrialSelect> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialSelectSiteOfType(std::istream& istr);
  virtual ~TrialSelectSiteOfType() {}

 protected:
  void serialize_trial_select_segment_(std::ostream& ostr) const;

 private:
  int site_type_;
};

inline std::shared_ptr<TrialSelectSiteOfType> MakeTrialSelectSiteOfType(
    const argtype &args = argtype()) {
  return std::make_shared<TrialSelectSiteOfType>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_SELECT_SITE_OF_TYPE_H_
