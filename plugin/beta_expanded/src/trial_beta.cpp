#include "beta_expanded/include/select_nothing.h"
#include "beta_expanded/include/perturb_beta.h"
#include "beta_expanded/include/compute_beta.h"
#include "beta_expanded/include/trial_beta.h"

namespace feasst {

std::shared_ptr<Trial> MakeTrialBeta(const argtype &args) {
  auto trial = std::make_shared<Trial>(args);
  trial->set_description("TrialBeta");
  trial->add_stage(
    MakeSelectNothing(args),
    MakePerturbBeta(args),
    args);
  trial->set(MakeComputeBeta());
  return trial;
}

}  // namespace feasst
