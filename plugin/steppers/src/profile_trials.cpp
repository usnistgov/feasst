#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
//#include "utils/include/timer.h"
#include "monte_carlo/include/trial_factory.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/profile_trials.h"

namespace feasst {

FEASST_MAPPER(ProfileTrials, argtype({{"trials_per_update", "1e3"}}));

ProfileTrials::ProfileTrials(argtype * args) : Analyze(args) {
  if (trials_per_update() < 1e3) {
    WARN("trials_per_update(" << trials_per_update() << ") should be " <<
      ">= 1e3 to ensure profiling does not slow the simulation.");
  }
}
ProfileTrials::ProfileTrials(argtype args) : ProfileTrials(&args) {
  feasst_check_all_used(args);
}

void ProfileTrials::initialize(MonteCarlo * mc) {
  Analyze::initialize(mc);
  profile_.clear();
  profile_.resize(mc->trial_factory().num());
  printer(header(*mc), output_file(mc->criteria()));
}

std::string ProfileTrials::header(const MonteCarlo& mc) const {
  std::stringstream ss;
  for (const std::shared_ptr<Trial>& trial : mc.trial_factory().trials()) {
    std::string name = trial->class_name();
    if (name == "Trial") {
      name = trial->description();
    }
    ss << name << ",";
  }
  ss << std::endl;
  return ss.str();
}

/// Store previous cpu time, then trigger update on the next step to compute
/// elapsed time.
void ProfileTrials::update(const MonteCarlo& mc) {
  DEBUG("is_previous " << is_previous_);
  if (is_previous_) {
    is_previous_ = false;
    const clock_t elapsed = clock() - previous_clock_;
    //const double elapsed = static_cast<double>(clock() - previous_clock_);
    DEBUG("elapsed " << elapsed);
    profile_[mc.trial_factory().last_index()].accumulate(elapsed);
  } else {
    is_previous_ = true;
    previous_clock_ = clock();
    trials_since_update_ = trials_per_update();
  }
}

std::string ProfileTrials::write(const MonteCarlo& mc) {
  std::stringstream ss;
  if (rewrite_header()) {
    ss << header(mc);
  }
  double total = 0.;
  for (const Accumulator& prof : profile_) {
    total += prof.sum();
  }
  for (const Accumulator& prof : profile_) {
    ss << prof.sum()/total << ",";
  }
  ss << std::endl;
  DEBUG(ss.str());
  return ss.str();
}

void ProfileTrials::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(6414, ostr);
  feasst_serialize_fstobj(profile_, ostr);
}

ProfileTrials::ProfileTrials(std::istream& istr) : Analyze(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 6414, "mismatch version:" << version);
  feasst_deserialize_fstobj(&profile_, istr);
}

}  // namespace feasst
