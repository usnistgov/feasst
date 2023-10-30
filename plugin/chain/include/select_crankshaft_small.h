
#ifndef FEASST_CHAIN_SELECT_CRANKSHAFT_SMALL_H_
#define FEASST_CHAIN_SELECT_CRANKSHAFT_SMALL_H_

#include "monte_carlo/include/trial_select_particle.h"

namespace feasst {

/**
  Similar to SelectSegment, but optimized for small molecules using anchors.
  The sites which move (mobile sites) are rotated about an axis defined by
  two anchors.
 */
class SelectCrankshaftSmall : public TrialSelectParticle {
 public:
  //@{
  /** @name Arguments
    - site[i]: add the (i+1)-th mobile site, beginning with i=1,
      where "site" is the first mobile site.
      The "[i]" is to be substituted for an integer 1, 2, 3 ...
    - anchor_site0: an anchor site to define rotation axis.
    - anchor_site1: a second, different anchor site to define rotation axis.
   */
  explicit SelectCrankshaftSmall(argtype args = argtype());
  explicit SelectCrankshaftSmall(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  void precompute(System * system) override {
    TrialSelectParticle::precompute(system); }

  bool select(const Select& perturbed,
    System* system,
    Random * random) override;

  std::shared_ptr<TrialSelect> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit SelectCrankshaftSmall(std::istream& istr);
  virtual ~SelectCrankshaftSmall() {}

  //@}
 protected:
  void serialize_select_crankshaft_small_(std::ostream& ostr) const;

 private:
};

inline std::shared_ptr<SelectCrankshaftSmall> MakeSelectCrankshaftSmall(
    argtype args = argtype()) {
  return std::make_shared<SelectCrankshaftSmall>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_SELECT_CRANKSHAFT_SMALL_H_
