#include "monte_carlo/include/trial_factory.h"

namespace feasst {

TrialFactory::TrialFactory() { class_name_ = "TrialFactory"; }

class MapTrialFactory {
 public:
  MapTrialFactory() {
    auto obj = MakeTrialFactory();
    obj->deserialize_map()["TrialFactory"] = obj;
  }
};

static MapTrialFactory mapper_ = MapTrialFactory();

void TrialFactory::add(std::shared_ptr<Trial> trial) {
  trials_.push_back(trial);

  // update probability of selection
  std::vector<double> weights;
  for (std::shared_ptr<Trial> trial : trials_) {
    weights.push_back(trial->weight());
  }
  cumulative_probability_ = cumulative_probability(weights);
  //std::stringstream ss;
  //ss << trials_.back()->class_name()"trial" << num_trials() - 1;
  // timer_.add(trials_.back()->class_name());
}

int TrialFactory::random_index(Random * random) {
  ASSERT(num_trials() > 0, "no trials to select");
  return random->index_from_cumulative_probability(cumulative_probability_);
}

bool TrialFactory::attempt(
    Criteria* criteria,
    System * system,
    int trial_index,
    Random * random) {
  increment_num_attempts();
  // timer_.start(0);
  if (num_trials() == 0) return false;
  if (trial_index == -1) {
    trial_index = random_index(random);
  }
  //timer_.start(index + 1);  // +1 for "other"
  const bool accepted = trials_[trial_index]->attempt(criteria, system, random);
  //timer_.end();
  if (accepted) increment_num_success_();
  return accepted;
}

void TrialFactory::revert(const int index,
                          const bool accepted,
                          System * system) {
  trials_[index]->revert(index, accepted, system);
  if (accepted) {
    decrement_num_success_();
  }
  decrement_num_attempts_();
}

void TrialFactory::mimic_trial_rejection(const int index) {
  DEBUG("index " << index << " " << trials_.size());
  trials_[index]->increment_num_attempts();
  increment_num_attempts();
}

std::string TrialFactory::status_header() const {
  std::stringstream ss;
  ss << "attempt ";
  for (int trial = 0; trial < num_trials(); ++trial) {
    ss << trials_[trial]->status_header() << " ";
  }
  return ss.str();
}

void TrialFactory::delay_finalize() {
  for (int trial = 0; trial < num_trials(); ++trial) {
    trials_[trial]->set_finalize_delayed(true);
  }
}

std::string TrialFactory::status() const {
  std::stringstream ss;
  ss << num_attempts() << " ";
  for (int trial = 0; trial < num_trials(); ++trial) {
    ss << trials_[trial]->status() << " ";
  }
  return ss.str();
}

void TrialFactory::reset_stats() {
  Trial::reset_stats();
  for (int trial = 0; trial < num_trials(); ++trial) {
    trials_[trial]->reset_stats();
  }
}

void TrialFactory::tune() {
  for (int trial = 0; trial < num_trials(); ++trial) {
    trials_[trial]->tune();
  }
}

void TrialFactory::precompute(Criteria * criteria, System * system) {
  for (std::shared_ptr<Trial> trial : trials_) {
    trial->precompute(criteria, system);
  }
}

bool TrialFactory::is_equal(const TrialFactory& factory) const {
  if (num_trials() != factory.num_trials()) {
    DEBUG("unequal number of trials: " << num_trials() << " "
      << factory.num_trials());
    return false;
  }
  if (!Trial::is_equal(factory)) {
    return false;
  }
  for (int it = 0; it < num_trials(); ++it) {
    if (!trials_[it]->is_equal(*(factory.trials_[it]))) {
      DEBUG("unequal trial" << it);
      return false;
    }
  }
  return true;
}

std::shared_ptr<Trial> TrialFactory::create(std::istream& istr) const {
  return std::make_shared<TrialFactory>(istr);
}

TrialFactory::TrialFactory(std::istream& istr) : Trial(istr) {
  // ASSERT(class_name_ == "TrialFactory", "name: " << class_name_);
  const int version = feasst_deserialize_version(istr);
  ASSERT(189 == version, "mismatch version: " << version);
  // HWH for unknown reasons, this function template does not work.
  // feasst_deserialize(&trials_, istr);
  { int dim1;
    istr >> dim1;
    trials_.resize(dim1);
    for (int index = 0; index < dim1; ++index) {
      // feasst_deserialize_fstobj((trials_)[index], istr);
      int existing;
      istr >> existing;
      if (existing != 0) {
        trials_[index] = trials_[index]->deserialize(istr);
      }
    }
  }
  feasst_deserialize(&cumulative_probability_, istr);
}

void TrialFactory::serialize_trial_factory_(std::ostream& ostr) const {
  serialize_trial_(ostr);
  feasst_serialize_version(189, ostr);
  feasst_serialize(trials_, ostr);
  feasst_serialize(cumulative_probability_, ostr);
}

void TrialFactory::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_trial_factory_(ostr);
}

}  // namespace feasst
