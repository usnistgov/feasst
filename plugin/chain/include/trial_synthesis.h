
#ifndef FEASST_CHAIN_TRIAL_SYNTHESIS_H_
#define FEASST_CHAIN_TRIAL_SYNTHESIS_H_

#include <string>
#include <memory>
#include "monte_carlo/include/trial.h"
#include "monte_carlo/include/trial_select_particle.h"
#include "monte_carlo/include/perturb_add.h"
#include "chain/include/trial_select_site_of_type.h"
#include "chain/include/perturb_site_type.h"
#include "chain/include/compute_synthesis.h"

namespace feasst {

/**
  Perform a synthesis reaction: R -> R' + P where
  a site type on R is changed and particle P is created.
 */
class TrialSynthesis : public Trial {
 public:
  /**
    args:
    - reactant_type: type of particle which is reacting, R.
    - reactant_site_type: type of site which is reacting on R.
    - new_site_type: type of site to create on R upon reaction completition.
    - product_type: type of particle, P, to add upon reaction.
   */
  TrialSynthesis(const argtype& args = argtype()) : Trial(args) {
    class_name_ = "TrialSynthesis";
    set(MakeComputeSynthesis());
    Arguments args_(args);
    args_.dont_check();
    const int reactant_type = args_.key("reactant_type").integer();
    const int reactant_site_type = args_.key("reactant_site_type").integer();
    const int new_site_type = args_.key("new_site_type").integer();
    const int product_type = args_.key("product_type").integer();
    DEBUG("product_type " << product_type);
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
        {"particle_type", str(product_type)}
      }),
      MakePerturbAdd()
    );
  }

  std::shared_ptr<Trial> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit TrialSynthesis(std::istream& istr);
  virtual ~TrialSynthesis() {}
};

inline std::shared_ptr<TrialSynthesis> MakeTrialSynthesis(
    const argtype &args = argtype()) {
  return std::make_shared<TrialSynthesis>(args);
}

}  // namespace feasst

#endif  // FEASST_CHAIN_TRIAL_SYNTHESIS_H_
