
#ifndef FEASST_CHAIN_TRIAL_SWAP_SITES_H_
#define FEASST_CHAIN_TRIAL_SWAP_SITES_H_

#include <string>
#include <memory>
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_compute_move.h"
#include "chain/include/trial_select_site_of_type.h"
#include "chain/include/perturb_site_type.h"

namespace feasst {

/**
  Swap the types of two sites in a particle.
 */
class TrialSwapSites : public Trial {
 public:
  /**
    args:
    - site_type1: type of site to swap.
    - site_type2: type of other site to swap.
   */
  TrialSwapSites(const argtype& args = argtype()) : Trial(args) {
    class_name_ = "TrialSwapSites";
    set(MakeTrialComputeMove());
    Arguments args_(args);
    args_.dont_check();
    const int site_type1 = args_.key("site_type1").integer();
    const int site_type2 = args_.key("site_type2").integer();
    ASSERT(site_type1 != site_type2, "site types should not match: " <<
      site_type1 << " " << site_type2);
    add_stage(
      MakeTrialSelectSiteOfType({{"site_type", str(site_type1)}}),
      MakePerturbSiteType({{"type", str(site_type2)}})
    );
    add_stage(
      MakeTrialSelectSiteOfType({{"site_type", str(site_type2)}}),
      MakePerturbSiteType({{"type", str(site_type1)}})
    );
  }

  std::shared_ptr<Trial> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialSwapSites(std::istream& istr);
  virtual ~TrialSwapSites() {}
};

inline std::shared_ptr<TrialSwapSites> MakeTrialSwapSites(
    const argtype &args = argtype()) {
  return std::make_shared<TrialSwapSites>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_SWAP_SITES_H_
