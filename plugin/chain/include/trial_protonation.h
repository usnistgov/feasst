
#ifndef FEASST_CHAIN_TRIAL_PROTONATION_H_
#define FEASST_CHAIN_TRIAL_PROTONATION_H_

#include <string>
#include <memory>
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/perturb_remove.h"
#include "chain/include/trial_select_site_of_type.h"
#include "chain/include/perturb_site_type.h"
#include "chain/include/compute_protonation.h"

namespace feasst {

/**
  Perform a protonation reaction: C(-) + P(+) -> CH where
  a negatively charged site on C reactions with a positively charged particle P
  P can be a proton or counter charge from salt via combination of
  deprotonation and grand canonical insertion of salt.
  See: https://doi.org/10.1063/1.4757284
 */
class TrialProtonation : public Trial {
 public:
  /**
    args:
    - reactant_type: type of particle which is reacting, C.
    - reactant_site_type: type of site which is reacting on C.
    - new_site_type: type of site to create on C upon reaction completition.
    - remove_type: type of particle, P, to remove upon reaction.
   */
  TrialProtonation(const argtype& args = argtype()) : Trial(args) {
    class_name_ = "TrialProtonation";
    set(MakeComputeProtonation());
    Arguments args_(args);
    args_.dont_check();
    const int reactant_type = args_.key("reactant_type").integer();
    const int reactant_site_type = args_.key("reactant_site_type").integer();
    const int new_site_type = args_.key("new_site_type").integer();
    const int remove_type = args_.key("remove_type").integer();
    DEBUG("remove_type " << remove_type);
    ASSERT(new_site_type != reactant_site_type, "site types should not match: " <<
      new_site_type << " " << reactant_site_type);
    DEBUG("adding first stage");
    add_stage(
      MakeTrialSelectSiteOfType({
        {"particle_type", str(reactant_type)},
        {"site_type", str(reactant_site_type)}
      }),
      MakePerturbSiteType({{"type", str(new_site_type)}})
    );
    DEBUG("adding second stage");
    add_stage(
      MakeTrialSelectParticle({
        {"particle_type", str(remove_type)}
      }),
      MakePerturbRemove()
    );
  }

  std::shared_ptr<Trial> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialProtonation(std::istream& istr);
  virtual ~TrialProtonation() {}
};

inline std::shared_ptr<TrialProtonation> MakeTrialProtonation(
    const argtype &args = argtype()) {
  return std::make_shared<TrialProtonation>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_PROTONATION_H_
