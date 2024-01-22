#include "utils/include/serialize.h"
#include "math/include/random.h"
#include "model_expanded/include/perturb_model.h"

namespace feasst {

PerturbModel::PerturbModel(argtype args) : PerturbModel(&args) {
  FEASST_CHECK_ALL_USED(args);
}
PerturbModel::PerturbModel(argtype * args) : Perturb(args) {
  class_name_ = "PerturbModel";
  potential_index_ = integer("potential_index", args, 0);
  disable_tunable_();
}

class MapPerturbModel {
 public:
  MapPerturbModel() {
    auto obj = MakePerturbModel();
    obj->deserialize_map()["PerturbModel"] = obj;
  }
};

static MapPerturbModel mapper_ = MapPerturbModel();

std::shared_ptr<Perturb> PerturbModel::create(std::istream& istr) const {
  return std::make_shared<PerturbModel>(istr);
}

void PerturbModel::change_model_index(const int delta_model, System * system) {
  Potential * pot = system->get_potential(potential_index_);
  pot->set_model_index(pot->model().model_index() + delta_model);
}

void PerturbModel::perturb(
    System * system,
    TrialSelect * select,
    Random * random,
    const bool is_position_held) {
  previous_model_ = system->potential(potential_index_).model().model_index();
  int delta_model = 1;
  if (random->coin_flip()) {
    delta_model = -1;
  }
  change_model_index(delta_model, system);
  DEBUG("changing model: " << previous_model_ << " by " << delta_model);
  set_revert_possible(true, select);
  select->set_trial_state(1);
}

void PerturbModel::revert(System * system) {
  DEBUG("revert_possible " << revert_possible());
  if (revert_possible()) {
    DEBUG("reverting model back to " << previous_model_);
    Potential * potential = system->get_potential(potential_index_);
    const int model_index = potential->model().model_index();
    change_model_index(previous_model_ - model_index, system);
  }
}

PerturbModel::PerturbModel(std::istream& istr) : Perturb(istr) {
  ASSERT(class_name_ == "PerturbModel", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(8935 == version, "mismatch version: " << version);
  feasst_deserialize(&potential_index_, istr);
}

void PerturbModel::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_perturb_(ostr);
  feasst_serialize_version(8935, ostr);
  feasst_serialize(potential_index_, ostr);
}

}  // namespace feasst
