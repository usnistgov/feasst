#include <string>
#include <memory>
#include "utils/include/serialize.h"
#include "monte_carlo/include/trial.h"
#include "math/include/random.h"

namespace feasst {

Trial::Trial(const argtype& args) {
  Arguments args_(args);
  args_.dont_check();
  set_new_only();
  set_finalize_delayed();
  weight_ = args_.key("weight").dflt("1").dble();
}

void Trial::add_stage(
  std::shared_ptr<TrialSelect> select,
  std::shared_ptr<Perturb> perturb,
  const argtype& args) {
  auto stage = std::make_shared<TrialStage>(args);
  stage->set(select);
  stage->set(perturb);
  add_(stage);
}

void Trial::set(const int index, std::shared_ptr<TrialStage> stage) {
  stages_[index] = stage;
  stages_ptr_[index] = stage.get();
}

void Trial::reset_stats() {
  DEBUG("reset_stats");
  num_attempts_ = 0;
  num_success_ = 0;
}

std::string Trial::status_header() const {
  std::stringstream ss;
  ss << "," << class_name_;
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
  if (num_attempts_ > 0) {
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
}

void Trial::revert(const int index, const bool accepted, System * system) {
  if (accepted) {
    revert(system);
    --num_success_;
  }
  --num_attempts_;
}

void Trial::finalize(System * system) {
  for (int index = num_stages() - 1; index >= 0; --index) {
    stages_[index]->finalize(system);
  }
}

bool Trial::attempt(Criteria * criteria, System * system, Random * random) {
  DEBUG("**********************************************************");
  DEBUG("* " << class_name() << " attempt " << num_attempts_ << " *");
  DEBUG("**********************************************************");
  ++num_attempts_;
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
  if (!acceptance_.reject()) {
    for (TrialStage * stage : stages_ptr_) {
      stage->set_mobile_physical(false, system);
    }
    compute_->perturb_and_acceptance(
      criteria, system, &acceptance_, &stages_ptr_, random);
  }
  DEBUG("num attempts: " << num_attempts_);
  if (criteria->is_accepted(acceptance_, *system, random->uniform())) {
    DEBUG("accepted");
    ++num_success_;
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

void Trial::serialize(std::ostream& ostr) const { FATAL("not implemented"); }

std::shared_ptr<Trial> Trial::create(std::istream& istr) const {
  FATAL("not implemented");
}

std::shared_ptr<Trial> Trial::deserialize(std::istream& istr) {
  return template_deserialize(deserialize_map(), istr,
    // true argument denotes rewinding to reread class name
    // this allows derived class constructor to read class name.
    true);
}

void Trial::refresh_stages_ptr_() {
  const int num = static_cast<int>(stages_.size());
  stages_ptr_.resize(num);
  for (int i = 0; i < num; ++i) {
    stages_ptr_[i] = stages_[i].get();
  }
}

bool Trial::is_equal(const Trial& trial) const {
  if (num_attempts_ != trial.num_attempts_) {
    INFO("unequal number of attempts:" << num_attempts_ << " "
      << trial.num_attempts_);
    return false;
  }
  if (num_success_ != trial.num_success_) {
    INFO("unequal number of success:" << num_success_ << " "
      << trial.num_success_);
    return false;
  }
  if (weight_ != trial.weight_) {
    INFO("unequal weight:" << weight_ << " " << trial.weight_);
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
  feasst_serialize(num_attempts_, ostr);
  feasst_serialize(num_success_, ostr);
  feasst_serialize(is_finalize_delayed_, ostr);
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
  feasst_deserialize(&num_attempts_, istr);
  feasst_deserialize(&num_success_, istr);
  feasst_deserialize(&is_finalize_delayed_, istr);
}

const std::vector<std::shared_ptr<Trial> >& Trial::trials() const {
  FATAL("not implemented");
}

const Trial& Trial::trial(const int index) const {
  FATAL("not implemented");
}

}  // namespace feasst
