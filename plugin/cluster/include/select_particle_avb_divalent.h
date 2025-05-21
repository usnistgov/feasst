
#ifndef FEASST_PARTICLE_AVB_SELECT_PARTICLE_AVB_DIVALENT_H_
#define FEASST_PARTICLE_AVB_SELECT_PARTICLE_AVB_DIVALENT_H_

#include <vector>
#include <memory>
#include "configuration/include/select.h"
#include "monte_carlo/include/trial_select_particle.h"

namespace feasst {

typedef std::map<std::string, std::string> argtype;

/**
  Using the first particle in perturbed, select another particle that has
  a site inside the aggregation volume (AV) of the target site.
 */
class SelectParticleAVBDivalent : public TrialSelect {
 public:
  //@{
  /** @name Arguments
   */

  /**
    args:
    - target_site_name: name of target site (default: 0).
    - site_name: name of site on particle_type to put in AV of target site
      (default: 0).
    - ghost: true if selecting a ghost (for adding).
    - neighbor_index: NeighborCriteria index contained in System (default: 0).
    - TrialSelect arguments.
   */
  explicit SelectParticleAVBDivalent(argtype args = argtype());

  //@}
  /** @name Public Functions
   */
  //@{

  void precompute(System * system) override;
  bool select(const Select& perturbed,
              System* system,
              Random * random,
              TrialSelect * previous_select) override;

  std::shared_ptr<TrialSelect> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit SelectParticleAVBDivalent(std::istream& istr);
  virtual ~SelectParticleAVBDivalent() {}

 protected:
  void serialize_select_particle_avb_divalent_(std::ostream& ostr) const;

  //@}
 private:
  int neighbor_;
  TrialSelectParticle select_mobile_;
  std::string target_site_name_;

  // temporary
  Select neighbors_;
};

inline std::shared_ptr<SelectParticleAVBDivalent> MakeSelectParticleAVBDivalent(
    const argtype &args = argtype()) {
  return std::make_shared<SelectParticleAVBDivalent>(args);
}

}  // namespace feasst

#endif  // FEASST_PARTICLE_AVB_SELECT_PARTICLE_AVB_DIVALENT_H_
