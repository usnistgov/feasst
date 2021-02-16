
#ifndef FEASST_PARTICLE_AVB_SELECT_PARTICLE_AVB_H_
#define FEASST_PARTICLE_AVB_SELECT_PARTICLE_AVB_H_

#include <vector>
#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial_select_particle.h"

namespace feasst {

/**
  Select a random particle as target, then selection another particle that has
  index of site either inside or outside the aggregation volume (AV) of
  target_site.
  Inside selections may optionally set the anchor to a second target particle,
  allowing for inside->inside moves.
 */
class SelectParticleAVB : public TrialSelect {
 public:
  /**
    args:
    - target_particle_type: type of target particle (default: 0).
    - target_site: index of target site (default: 0).
    - site: index of site on particle_type to put in AV of target site
      (default: 0).
    - grand_canonical: true if used for grand canonical, false otherwise.
    - inside: true if selecting in the AV, otherwise out (default: true).
      Not implemented for grand_canonical.
    - second_target: if true, set anchor to a second target particle
      (default: false).
    - neighbor_index: NeighborCriteria index contained in System (default: 0).
   */
  SelectParticleAVB(argtype args = argtype());
  SelectParticleAVB(argtype * args);

  // ghost is set in precompute
  void precompute(System * system) override;
  bool select(const Select& perturbed,
              System* system,
              Random * random) override;

  std::shared_ptr<TrialSelect> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit SelectParticleAVB(std::istream& istr);
  virtual ~SelectParticleAVB() {}

 protected:
  void serialize_select_particle_avb_(std::ostream& ostr) const;

 private:
  int neighbor_;
  int site_index_;
  bool grand_canonical_;
  bool inside_;
  bool is_second_target_;
  TrialSelectParticle select_target_;
  TrialSelectParticle select_mobile_;

  // temporary
  Select target_;
  Select second_target_;
  Select empty_;
  Select neighbors_;
};

inline std::shared_ptr<SelectParticleAVB> MakeSelectParticleAVB(
    argtype args = argtype()) {
  return std::make_shared<SelectParticleAVB>(args);
}

}  // namespace feasst

#endif  // FEASST_PARTICLE_AVB_SELECT_PARTICLE_AVB_H_
