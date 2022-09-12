#include "utils/include/serialize.h"
//#include "utils/include/timer.h"
#include "steppers/include/profile_trials.h"

namespace feasst {

class MapProfileTrials {
 public:
  MapProfileTrials() {
    auto obj = MakeProfileTrials();
    obj->deserialize_map()["ProfileTrials"] = obj;
  }
};

static MapProfileTrials mapper_ = MapProfileTrials();

ProfileTrials::ProfileTrials(argtype * args) : Analyze(args) {}
ProfileTrials::ProfileTrials(argtype args) : ProfileTrials(&args) {
  FEASST_CHECK_ALL_USED(args);
}

void ProfileTrials::initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  profile_.clear();
  profile_.resize(trial_factory->num());
  printer(header(*criteria, *system, *trial_factory),
          file_name(*criteria));
}

std::string ProfileTrials::header(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) const {
  std::stringstream ss;
  for (const std::shared_ptr<Trial>& trial : trial_factory.trials()) {
    ss << trial->class_name() << " ";
  }
  ss << std::endl;
  return ss.str();
}

/// Store previous cpu time, then trigger update on the next step to compute
/// elapsed time.
void ProfileTrials::update(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  DEBUG("is_previous " << is_previous_);
  if (is_previous_) {
    is_previous_ = false;
    const clock_t elapsed = clock() - previous_clock_;
    //const double elapsed = static_cast<double>(clock() - previous_clock_);
    DEBUG("elapsed " << elapsed);
    profile_[trial_factory.last_index()].accumulate(elapsed);
  } else {
    is_previous_ = true;
    previous_clock_ = clock();
    trials_since_update_ = trials_per_update();
  }
}

std::string ProfileTrials::write(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  std::stringstream ss;
  for (const Accumulator& prof : profile_) {
    ss << prof.sum() << " ";
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
