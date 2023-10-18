#ifndef FEASST_CHAIN_TRIAL_SWAP_SWITES_H_
#define FEASST_CHAIN_TRIAL_SWAP_SWITES_H_

#include <memory>
#include "utils/include/arguments.h"
#include "monte_carlo/include/trial.h"

namespace feasst {

/**
  Swap the types of two sites in a particle.
 */
class TrialSwapSites : public Trial {
 public:
  //@{
  /** @name Arguments
   */

  /**
    args:
    - site_type1: type of site to swap.
    - site_type2: type of other site to swap.
    - Trial arguments.
    - TrialStage arguments.
    - TrialSelect arguments.
   */
  explicit TrialSwapSites(argtype args = argtype());
  explicit TrialSwapSites(argtype * args);

  //@}
  /** @name Public Functions
   */
  //@{

  std::shared_ptr<Trial> create(std::istream& istr) const override {
    return std::make_shared<TrialSwapSites>(istr); }
  std::shared_ptr<Trial> create(argtype * args) const override {
    return std::make_shared<TrialSwapSites>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit TrialSwapSites(std::istream& istr);
  virtual ~TrialSwapSites() {}
  //@}
};

inline std::shared_ptr<TrialSwapSites> MakeTrialSwapSites(argtype args = argtype()) {
  return std::make_shared<TrialSwapSites>(args); }

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_SWAP_SWITES_H_
