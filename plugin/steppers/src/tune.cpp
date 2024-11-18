#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "monte_carlo/include/trial_factory.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/tunable.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/perturb.h"
#include "monte_carlo/include/acceptance.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/tune.h"

namespace feasst {

FEASST_MAPPER(Tune,);

Tune::Tune(argtype * args) : Modify(args) {
  trials_per_tune_ = integer("trials_per_tune", args, 1e3);
  ASSERT(trials_per_update() == 1, "requires 1 trial per update");
  if (trials_per_write() == 1) {
    set_trials_per_write(-1);
  }
}
Tune::Tune(argtype args) : Tune(&args) { feasst_check_all_used(args); }

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

void Tune::initialize(MonteCarlo * mc) {
  const TrialFactory& tfac = mc->trial_factory();
  const int num_trials = tfac.num();
  values_.resize(num_trials);
  num_attempts_.resize(num_trials, 0);
  num_accepted_.resize(num_trials, 0);
  for (int trial = 0; trial < num_trials; ++trial) {
    if (tfac.trial(trial).num_stages() > 0) {
      values_[trial] = tfac.trial(trial).stage(0).perturb().tunable().value();
    }
  }
  if (trials_per_write() != -1 && !output_file().empty()) {
    printer(header(*mc), output_file(mc->criteria()));
  }
}

void Tune::update(MonteCarlo * mc) {
  TrialFactory * tfac = mc->get_trial_factory();
  const int trial = tfac->last_index();
  DEBUG("last trial attempted: " << trial);
  if (trial < min_num(*tfac)) {
    if (!tfac->trial(trial).accept().reject()) {
      // update acceptance statistics
      int * num_attempts = &num_attempts_[trial];
      *num_attempts += 1;
      DEBUG("num_attempts: " << *num_attempts);
      int * num_accepted = &num_accepted_[trial];
      if (mc->criteria().was_accepted()) {
        *num_accepted += 1;
      }

      // check for tuning
      if (tfac->trial(trial).num_stages() > 0) {
        const Tunable& tunable =
          tfac->trial(trial).stage(0).perturb().tunable();
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
          tfac->set_tunable(trial, *value);
        }
      }
    }
  }
}

std::string Tune::header(const MonteCarlo& mc) const {
  const TrialFactory& tfac = mc.trial_factory();
  std::stringstream ss;
  for (int trial = 0; trial < min_num(tfac); ++trial) {
    ss << tfac.trial(trial).class_name() << ",acceptance,";
  }
  ss << std::endl;
  return ss.str();
}

std::string Tune::write(MonteCarlo * mc) {
  std::stringstream ss;
  if (rewrite_header()) {
    ss << header(*mc);
  }
  for (int trial = 0; trial < min_num(mc->trial_factory()); ++trial) {
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
