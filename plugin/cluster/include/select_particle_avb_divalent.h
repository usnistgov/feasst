
#ifndef FEASST_PARTICLE_AVB_SELECT_PARTICLE_AVB_DIVALENT_H_
#define FEASST_PARTICLE_AVB_SELECT_PARTICLE_AVB_DIVALENT_H_

#include <vector>
#include <memory>
#include "utils/include/arguments.h"
#include "system/include/neighbor_criteria.h"
#include "monte_carlo/include/trial_select_particle.h"

namespace feasst {

/**
  Using the first particle in perturbed, select another particle that has
  site_index inside the aggregation volume (AV) of target_site_index.
 */
class SelectParticleAVBDivalent : public TrialSelect {
 public:
  /**
    args:
    - target_site_index: index of target site (default: 0).
    - site_index: index of site on particle_type to put in AV of target site
      (default: 0).
    - ghost: true if selecting a ghost (for adding).
   */
  SelectParticleAVBDivalent(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype& args = argtype());

  bool select(const Select& perturbed,
              System* system,
              Random * random) override;

  std::shared_ptr<TrialSelect> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit SelectParticleAVBDivalent(std::istream& istr);
  virtual ~SelectParticleAVBDivalent() {}

 protected:
  void serialize_select_particle_avb_divalent_(std::ostream& ostr) const;

 private:
  std::shared_ptr<NeighborCriteria> neighbor_criteria_;
  TrialSelectParticle select_mobile_;

  // temporary
  Select neighbors_;
};

inline std::shared_ptr<SelectParticleAVBDivalent> MakeSelectParticleAVBDivalent(
    std::shared_ptr<NeighborCriteria> neighbor_criteria,
    const argtype &args = argtype()) {
  return std::make_shared<SelectParticleAVBDivalent>(neighbor_criteria, args);
}

}  // namespace feasst

#endif  // FEASST_PARTICLE_AVB_SELECT_PARTICLE_AVB_DIVALENT_H_
