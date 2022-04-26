#include "steppers/include/tune_per_state.h"
#include "utils/include/serialize.h"

namespace feasst {

class MapTunePerState {
 public:
  MapTunePerState() {
    TunePerState().deserialize_map()["TunePerState"] = MakeTunePerState();
  }
};

static MapTunePerState mapper_ = MapTunePerState();

TunePerState::TunePerState(argtype * args) : Modify(args) {
  trials_per_tune_ = integer("trials_per_tune", args, 1e2);
  ASSERT(trials_per_update() == 1, "requires 1 trial per update");
  stop_after_iteration_ = integer("stop_after_iteration", args, -1);
}
TunePerState::TunePerState(argtype args) : TunePerState(&args) { check_all_used(args); }

void TunePerState::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(256, ostr);
  feasst_serialize(trials_per_tune_, ostr);
  feasst_serialize(stop_after_iteration_, ostr);
  feasst_serialize(values_, ostr);
  feasst_serialize(num_attempts_, ostr);
  feasst_serialize(num_accepted_, ostr);
}

TunePerState::TunePerState(std::istream& istr) : Modify(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 256, "version mismatch:" << version);
  feasst_deserialize(&trials_per_tune_, istr);
  feasst_deserialize(&stop_after_iteration_, istr);
  feasst_deserialize(&values_, istr);
  feasst_deserialize(&num_attempts_, istr);
  feasst_deserialize(&num_accepted_, istr);
}

void TunePerState::initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  const int num_trials = trial_factory->num();
  values_.resize(num_trials);
  num_attempts_.resize(num_trials);
  num_accepted_.resize(num_trials);
  for (int trial = 0; trial < num_trials; ++trial) {
    const double initial_value =
      trial_factory->trial(trial).stage(0).perturb().tunable().value();
    DEBUG("initial_value: " << initial_value);
    values_[trial].resize(criteria->num_states(), initial_value);
    num_attempts_[trial].resize(criteria->num_states(), 0);
    num_accepted_[trial].resize(criteria->num_states(), 0);
  }
}

void TunePerState::update(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  //WARN("hwh also, adjusting windows.. carry over tune data? .. or not");
  for (int trial = 0; trial < static_cast<int>(values_.size()); ++trial) {
  //for (int trial = 0; trial < trial_factory->num(); ++trial) {
    const Tunable& tunable =
      trial_factory->trial(trial).stage(0).perturb().tunable();
    if (tunable.is_enabled()) {
      const int current_state = criteria->state();
      DEBUG("current_state: " << current_state);
      double * value = &values_[trial][current_state];

      // if trial is the last one attempted, update the acceptance stats
      // then check for tuning
      if (trial == trial_factory->last_index() &&
          (stop_after_iteration_ == -1 ||
           stop_after_iteration_ > criteria->num_iterations())) {
        int * num_attempts = &num_attempts_[trial][current_state];
        *num_attempts += 1;
        //++num_attempts_[trial][current_state];
        DEBUG("num_attempts: " << *num_attempts);
        int * num_accepted = &num_accepted_[trial][current_state];
        if (criteria->was_accepted()) {
          *num_accepted += 1;
        }
        DEBUG("num_accepted: " << *num_accepted);
        DEBUG("trials_per_tune_: " << trials_per_tune_);
        if (*num_attempts == trials_per_tune_) {
          const double acceptance = *num_accepted/
                static_cast<double>(*num_attempts);
          DEBUG("acceptance: " << acceptance);
          double val = *value;
          val *= 1 + tunable.percent_change()*(acceptance - tunable.target());
          if (!tunable.is_bound() ||
              (val <= tunable.max() && val >= tunable.min())) {
            *value = val;
          }
          *num_accepted = 0;
          *num_attempts = 0;
        }
      }
      trial_factory->set_tunable(trial, *value);
    }
  }
}

std::string TunePerState::write(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  std::stringstream ss;
  //for (int trial = 0; trial < trial_factory->num(); ++trial) {
  for (int trial = 0; trial < static_cast<int>(values_.size()); ++trial) {
    ss << trial_factory->trial(trial).class_name() << ",";
  }
  ss << std::endl;
  for (int state = 0; state < criteria->num_states(); ++state) {
    //for (int trial = 0; trial < trial_factory->num(); ++trial) {
    for (int trial = 0; trial < static_cast<int>(values_.size()); ++trial) {
      ss << values_[trial][state] << ",";
    }
    ss << std::endl;
  }
  return ss.str();
}

}  // namespace feasst
