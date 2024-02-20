#include "utils/include/serialize.h"
#include "math/include/random.h"
#include "beta_expanded/include/perturb_beta.h"

namespace feasst {

PerturbBeta::PerturbBeta(argtype args) : PerturbBeta(&args) {
  FEASST_CHECK_ALL_USED(args);
}
PerturbBeta::PerturbBeta(argtype * args) : Perturb(args) {
  class_name_ = "PerturbBeta";
  fixed_beta_change_ = dble("fixed_beta_change", args);
  disable_tunable_();
}

class MapPerturbBeta {
 public:
  MapPerturbBeta() {
    auto obj = MakePerturbBeta({{"fixed_beta_change", "1"}});
    obj->deserialize_map()["PerturbBeta"] = obj;
  }
};

static MapPerturbBeta mapper_ = MapPerturbBeta();

std::shared_ptr<Perturb> PerturbBeta::create(std::istream& istr) const {
  return std::make_shared<PerturbBeta>(istr);
}

void PerturbBeta::change_beta(const double delta_beta, System * system) {
  system->set_beta(system->thermo_params().beta() + delta_beta);
}

void PerturbBeta::perturb(
    System * system,
    TrialSelect * select,
    Random * random,
    const bool is_position_held,
    Acceptance * acceptance) {
  previous_beta_ = system->thermo_params().beta();
  double delta_beta = fixed_beta_change_;
  if (random->coin_flip()) delta_beta *= -1;
  change_beta(delta_beta, system);
  DEBUG("changing beta: " << previous_beta_ << " by " << delta_beta);
  set_revert_possible(true, select);
  select->set_trial_state(1);
}

void PerturbBeta::revert(System * system) {
  DEBUG("revert_possible " << revert_possible());
  if (revert_possible()) {
    DEBUG("reverting beta back to " << previous_beta_);
    change_beta(previous_beta_ - system->thermo_params().beta(), system);
  }
}

PerturbBeta::PerturbBeta(std::istream& istr)
  : Perturb(istr) {
  ASSERT(class_name_ == "PerturbBeta", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(1634 == version, "mismatch version: " << version);
  feasst_deserialize(&fixed_beta_change_, istr);
}

void PerturbBeta::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_(ostr);
  feasst_serialize_version(1634, ostr);
  feasst_serialize(fixed_beta_change_, ostr);
}

}  // namespace feasst
