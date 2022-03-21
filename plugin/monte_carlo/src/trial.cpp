#include <string>
#include <memory>
#include "utils/include/serialize.h"
#include "monte_carlo/include/trial.h"
#include "math/include/random.h"

namespace feasst {

class MapTrial {
 public:
  MapTrial() {
    auto obj = MakeTrial();
    obj->deserialize_map()["Trial"] = obj;
  }
};

static MapTrial mapper_ = MapTrial();

Trial::Trial(argtype args) : Trial(&args) { check_all_used(args); }
Trial::Trial(argtype * args) {
  set_finalize_delayed();
  weight_ = dble("weight", args, 1);
  data_.get_int64_1D()->resize(3);
  reset_stats();
}

void Trial::add_stage(
  std::shared_ptr<TrialSelect> select,
  std::shared_ptr<Perturb> perturb,
  argtype * stage_args) {
  auto stage = std::make_shared<TrialStage>(stage_args);
  stage->set(select);
  stage->set(perturb);
  add_(stage);
}

void Trial::add_stage(std::shared_ptr<TrialSelect> select,
                      std::shared_ptr<Perturb> perturb) {
  argtype empty_args;
  add_stage(select, perturb, &empty_args);
}

void Trial::set(const int index, std::shared_ptr<TrialStage> stage) {
  stages_[index] = stage;
  stages_ptr_[index] = stage.get();
}

void Trial::reset_stats() {
  DEBUG("reset_stats");
  *num_attempts_() = 0;
  *num_success_() = 0;
  *num_auto_reject_() = 0;
}

std::string Trial::status_header() const {
  std::stringstream ss;
  ss << ",";
  if (class_name_ != "Trial") {
    ss << class_name_;
  } else {
    ss << description();
  }
  for (const TrialStage * stage : stages_ptr_) {
    ss << stage->status_header();
  }
  return ss.str();
}

std::string Trial::status() const {
  std::stringstream ss;
  ss << "," << acceptance();
  for (const TrialStage * stage : stages_ptr_) {
    ss << stage->status();
  }
  return ss.str();
}

void Trial::tune() {
  int num_real_attempts = num_attempts() - num_auto_reject();
  DEBUG("num " << num_attempts());
  DEBUG("num_auto_rej " << num_auto_reject());
  DEBUG("num_real " << num_real_attempts);
  if (num_real_attempts > 0) {
    for (auto stage : stages_) stage->tune(acceptance());
    reset_stats();
  }
}

void Trial::precompute(Criteria * criteria, System * system) {
  for (std::shared_ptr<TrialStage> stage : stages_) {
    stage->precompute(system);
  }
}

void Trial::revert(System * system) {
  for (int index = num_stages() - 1; index >= 0; --index) {
    stages_[index]->revert(system);
  }
  system->revert(acceptance_.perturbed());
}

void Trial::revert(const int index,
    const bool accepted,
    const bool auto_rejected,
    System * system) {
  if (accepted) {
    revert(system);
    decrement_num_success_();
  }
  //ASSERT(!auto_rejected, "er");
  DEBUG("auto_rejected " << auto_rejected);
  if (auto_rejected) *num_auto_reject_() -= 1;
  decrement_num_attempts_();
}

void Trial::finalize(System * system) {
  for (int index = num_stages() - 1; index >= 0; --index) {
    stages_[index]->finalize(system);
  }
  DEBUG("finalize perturbed");
  system->finalize(acceptance_.perturbed());
  DEBUG("done finalizing perturbed");
}

bool Trial::attempt(Criteria * criteria, System * system, Random * random) {
  DEBUG("**********************************************************");
  DEBUG("* " << class_name() << " " << description() << " attempt " << num_attempts() << " *");
  DEBUG("**********************************************************");
  DEBUG("num particles: " << system->configuration().num_particles());
  DEBUG("num ghosts: " << system->configuration().particles().num() -
                         system->configuration().num_particles());
  DEBUG("existing: " << system->configuration().group_select(0).str());
  DEBUG("num of type 0: " << system->configuration().num_particles_of_type(0));
  DEBUG("current_energy: " << criteria->current_energy());
  DEBUG("all: " << system->configuration().selection_of_all().str());
  increment_num_attempts();
  acceptance_.reset();
  criteria->before_attempt(*system);
  before_select(&acceptance_, criteria);

  // Perform selections. If one selection fails, do not continue selecting.
  for (TrialStage * stage : stages_ptr_) {
    stage->before_select();
    DEBUG("selecting");
    if (!acceptance_.reject()) {
      stage->select(system, &acceptance_, random);
    }
  }
  if (acceptance_.reject()) {
    DEBUG("auto rejected at selection");
  } else {
    for (TrialStage * stage : stages_ptr_) {
      stage->set_mobile_physical(false, system);
    }
    compute_->perturb_and_acceptance(
      criteria, system, &acceptance_, &stages_ptr_, random);
  }
  DEBUG("num attempts: " << num_attempts());
  if (acceptance_.reject()) {
    DEBUG("auto reject");
    *num_auto_reject_() += 1;
  }
  if (criteria->is_accepted(*system, &acceptance_, random)) {
    DEBUG("accepted");
    increment_num_success_();
    if (!is_finalize_delayed_) {
      finalize(system);
    }
    return true;
  } else {
    DEBUG("rejected");
    revert(system);
    return false;
  }
}

std::map<std::string, std::shared_ptr<Trial> >& Trial::deserialize_map() {
  static std::map<std::string, std::shared_ptr<Trial> >* ans =
     new std::map<std::string, std::shared_ptr<Trial> >();
  return *ans;
}

void Trial::serialize(std::ostream& ostr) const {
  ostr << class_name() << " ";
  serialize_trial_(ostr);
}

std::shared_ptr<Trial> Trial::create(std::istream& istr) const {
  return std::make_shared<Trial>(istr);
}

std::shared_ptr<Trial> Trial::create(argtype * args) const {
  FATAL("not implemented");
}

std::shared_ptr<Trial> Trial::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

std::shared_ptr<Trial> Trial::factory(const std::string name, argtype * args) {
  return template_factory(deserialize_map(), name, args);
}

void Trial::refresh_stages_ptr_() {
  const int num = static_cast<int>(stages_.size());
  stages_ptr_.resize(num);
  for (int i = 0; i < num; ++i) {
    stages_ptr_[i] = stages_[i].get();
  }
}

bool Trial::is_equal(const Trial& trial) const {
  if (num_attempts() != trial.num_attempts()) {
    DEBUG("unequal number of attempts:" << num_attempts() << " "
      << trial.num_attempts());
    return false;
  }
  if (num_success() != trial.num_success()) {
    DEBUG("unequal number of success:" << num_success() << " "
      << trial.num_success());
    return false;
  }
  if (weight_ != trial.weight_) {
    DEBUG("unequal weight:" << weight_ << " " << trial.weight_);
    return false;
  }
  if (num_stages() > 0) {
    if (!stages_[0]->perturb().tunable().is_equal(
        trial.stages_[0]->perturb().tunable())) {
      return false;
    }
  }
  return true;
}

void Trial::serialize_trial_(std::ostream& ostr) const {
  feasst_serialize_version(570, ostr);
  feasst_serialize(stages_, ostr);
  // desererialize: refresh stages_ptr_
  feasst_serialize_fstdr(compute_, ostr);
  feasst_serialize(weight_, ostr);
  feasst_serialize(description_, ostr);
  //feasst_serialize(num_attempts_, ostr);
  //feasst_serialize(num_success_, ostr);
  feasst_serialize(is_finalize_delayed_, ostr);
  feasst_serialize_fstobj(data_, ostr);
}

Trial::Trial(std::istream& istr) {
  istr >> class_name_;
  const int version = feasst_deserialize_version(istr);
  ASSERT(570 == version, "mismatch version: " << version);
  // HWH for unknown reasons, this function template does not work.
  // feasst_deserialize(&stages_, istr);
  { int dim1;
    istr >> dim1;
    stages_.resize(dim1);
    for (int index = 0; index < dim1; ++index) {
      // feasst_deserialize_fstobj((stages_)[index], istr);
      int existing;
      istr >> existing;
      if (existing != 0) {
        stages_[index] = std::make_shared<TrialStage>(istr);
      }
    }
  }

  refresh_stages_ptr_();
  // HWH for unknown reasons, this function template does not work.
  //feasst_deserialize_fstdr(compute_, istr);
  { int existing;
    istr >> existing;
    if (existing != 0) {
      compute_ = compute_->deserialize(istr);
    }
  }
  feasst_deserialize(&weight_, istr);
  feasst_deserialize(&description_, istr);
  //feasst_deserialize(&num_attempts_, istr);
  //feasst_deserialize(&num_success_, istr);
  feasst_deserialize(&is_finalize_delayed_, istr);
  feasst_deserialize_fstobj(&data_, istr);
}

const std::vector<std::shared_ptr<Trial> >& Trial::trials() const {
  FATAL("not implemented");
}

const Trial& Trial::trial(const int index) const {
  FATAL("not implemented");
}

double Trial::acceptance() const {
  int num_real_attempts = num_attempts() - num_auto_reject();
  if (num_real_attempts == 0) return -1;
  return static_cast<double>(num_success())/
         static_cast<double>(num_real_attempts);
}

}  // namespace feasst
