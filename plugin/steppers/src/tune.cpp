#include "steppers/include/tune.h"
#include "utils/include/serialize.h"

namespace feasst {

class MapTune {
 public:
  MapTune() {
    Tune().deserialize_map()["Tune"] = MakeTune();
  }
};

static MapTune mapper_ = MapTune();

Tune::Tune(argtype * args) : Modify(args) {
  trials_per_tune_ = integer("trials_per_tune", args, 1e3);
  ASSERT(trials_per_update() == 1, "requires 1 trial per update");
}
Tune::Tune(argtype args) : Tune(&args) { FEASST_CHECK_ALL_USED(args); }

void Tune::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(256, ostr);
  feasst_serialize(trials_per_tune_, ostr);
  feasst_serialize(values_, ostr);
  feasst_serialize(num_attempts_, ostr);
  feasst_serialize(num_accepted_, ostr);
}

Tune::Tune(std::istream& istr) : Modify(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 256, "version mismatch:" << version);
  feasst_deserialize(&trials_per_tune_, istr);
  feasst_deserialize(&values_, istr);
  feasst_deserialize(&num_attempts_, istr);
  feasst_deserialize(&num_accepted_, istr);
}

void Tune::initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  const int num_trials = trial_factory->num();
  values_.resize(num_trials);
  num_attempts_.resize(num_trials, 0);
  num_accepted_.resize(num_trials, 0);
  for (int trial = 0; trial < num_trials; ++trial) {
    if (trial_factory->trial(trial).num_stages() > 0) {
      values_[trial] = trial_factory->trial(trial).stage(0).perturb().tunable().value();
    }
  }
  printer(header(*criteria, *system, *trial_factory),
          file_name(*criteria));
}

void Tune::update(Criteria * criteria,
    System * system,
    Random * random,
    TrialFactory * trial_factory) {
  const int trial = trial_factory->last_index();
  if (trial < min_num(*trial_factory)) {
    if (!trial_factory->trial(trial).accept().reject()) {
      // update acceptance statistics
      int * num_attempts = &num_attempts_[trial];
      *num_attempts += 1;
      DEBUG("num_attempts: " << *num_attempts);
      int * num_accepted = &num_accepted_[trial];
      if (criteria->was_accepted()) {
        *num_accepted += 1;
      }

      // check for tuning
      if (trial_factory->trial(trial).num_stages() > 0) {
        const Tunable& tunable =
          trial_factory->trial(trial).stage(0).perturb().tunable();
        if (tunable.is_enabled()) {
          double * value = &values_[trial];
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
          trial_factory->set_tunable(trial, *value);
        }
      }
    }
  }
}

std::string Tune::header(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) const {
  std::stringstream ss;
  for (int trial = 0; trial < min_num(trial_factory); ++trial) {
    ss << trial_factory.trial(trial).class_name() << ",acceptance,";
  }
  ss << std::endl;
  return ss.str();
}

std::string Tune::write(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  std::stringstream ss;
  if (rewrite_header()) {
    ss << header(*criteria, *system, *trial_factory);
  }
  for (int trial = 0; trial < min_num(*trial_factory); ++trial) {
    DEBUG("trial " << trial);
    DEBUG("valsz " << values_.size());
    DEBUG("attsz " << num_attempts_.size());
    DEBUG("accsz " << num_accepted_.size());
    DEBUG("val " << values_[trial]);
    DEBUG("att " << num_attempts_[trial]);
    DEBUG("acc " << num_accepted_[trial]);
    ss << values_[trial] << ",";
    if (num_attempts_[trial] == 0) {
      ss << "0,";
    } else {
      ss << num_accepted_[trial]/static_cast<double>(num_attempts_[trial]) << ",";
    }
  }
  ss << std::endl;
  return ss.str();
}

int Tune::min_num(const TrialFactory& trial_factory) const {
  return std::min(trial_factory.num(), static_cast<int>(values_.size()));
}

}  // namespace feasst
