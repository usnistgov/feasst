
#ifndef FEASST_CHAIN_TRIAL_DEPROTONATION_H_
#define FEASST_CHAIN_TRIAL_DEPROTONATION_H_

#include <string>
#include <memory>
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/perturb_add.h"
#include "chain/include/trial_select_site_of_type.h"
#include "chain/include/perturb_site_type.h"
#include "chain/include/compute_deprotonation.h"

namespace feasst {

/**
  Perform a deprotonation reaction: CH -> C(-) + P(+) where
  a site on C becomes negatively charged and a positively charged particle P
  is created.
  P can be a proton or counter charge from salt via combination of
  deprotonation and grand canonical insertion of salt.
  See: https://doi.org/10.1063/1.4757284
 */
class TrialDeprotonation : public Trial {
 public:
  /**
    args:
    - reactant_type: type of particle which is reacting, C.
    - reactant_site_type: type of site which is reacting on C.
    - new_site_type: type of site to create on R upon reaction completition.
    - add_type: type of particle, P, to add upon reaction.
   */
  TrialDeprotonation(const argtype& args = argtype()) : Trial(args) {
    class_name_ = "TrialDeprotonation";
    set(MakeComputeDeprotonation());
    Arguments args_(args);
    args_.dont_check();
    const int reactant_type = args_.key("reactant_type").integer();
    const int reactant_site_type = args_.key("reactant_site_type").integer();
    const int new_site_type = args_.key("new_site_type").integer();
    const int add_type = args_.key("add_type").integer();
    DEBUG("add_type " << add_type);
    ASSERT(new_site_type != reactant_site_type, "site types should not match: " <<
      new_site_type << " " << reactant_site_type);
    add_stage(
      MakeTrialSelectSiteOfType({
        {"particle_type", str(reactant_type)},
        {"site_type", str(reactant_site_type)}
      }),
      MakePerturbSiteType({{"type", str(new_site_type)}})
    );
    add_stage(
      MakeTrialSelectParticle({
        {"particle_type", str(add_type)}
      }),
      MakePerturbAdd()
    );
  }

  std::shared_ptr<Trial> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialDeprotonation(std::istream& istr);
  virtual ~TrialDeprotonation() {}
};

inline std::shared_ptr<TrialDeprotonation> MakeTrialDeprotonation(
    const argtype &args = argtype()) {
  return std::make_shared<TrialDeprotonation>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_DEPROTONATION_H_
