#include <cmath>
#include "utils/include/io.h"
#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "math/include/random.h"
#include "configuration/include/select.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/configuration.h"
#include "configuration/include/domain.h"
#include "system/include/system.h"
#include "system/include/thermo_params.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/acceptance.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/trial_compute_add.h"
#include "monte_carlo/include/trial_compute_remove.h"
#include "monte_carlo/include/trial_compute_add_remove.h"
#include "monte_carlo/include/trial_select.h"

namespace feasst {

TrialComputeAddRemove::TrialComputeAddRemove(argtype * args) : TrialCompute(args) {
  class_name_ = "TrialComputeAddRemove";
  argtype args2 = *args;
  add_ = std::make_unique<TrialComputeAdd>(args);
  rm_ = std::make_unique<TrialComputeRemove>(args);
}
TrialComputeAddRemove::TrialComputeAddRemove(argtype args) : TrialComputeAddRemove(&args) {
  feasst_check_all_used(args);
}
TrialComputeAddRemove::~TrialComputeAddRemove() {}

FEASST_MAPPER(TrialComputeAddRemove,);

void TrialComputeAddRemove::perturb_and_acceptance(
    Criteria * criteria,
    System * system,
    Acceptance * acceptance,
    std::vector<TrialStage*> * stages,
    Random * random) {
  DEBUG("TrialComputeAddRemove");
  TrialCompute * compute;
  if ((*stages)[0]->select().is_ghost()) {
    compute = add_.get();
  } else {
    compute = rm_.get();
  }
  compute->perturb_and_acceptance(criteria, system, acceptance, stages, random);
}

std::shared_ptr<TrialCompute> TrialComputeAddRemove::create(std::istream& istr) const {
  return std::make_shared<TrialComputeAddRemove>(istr);
}

TrialComputeAddRemove::TrialComputeAddRemove(std::istream& istr)
  : TrialCompute(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(5740 == version, "mismatch version: " << version);
  feasst_deserialize(add_, istr);
  feasst_deserialize(rm_, istr);
}

void TrialComputeAddRemove::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_compute_(ostr);
  feasst_serialize_version(5740, ostr);
  feasst_serialize(add_, ostr);
  feasst_serialize(rm_, ostr);
}

}  // namespace feasst
