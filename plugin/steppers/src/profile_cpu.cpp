#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "utils/include/timer_rdtsc.h"
#include "math/include/accumulator.h"
#include "monte_carlo/include/criteria.h"
#include "monte_carlo/include/analyze_factory.h"
#include "monte_carlo/include/modify_factory.h"
#include "monte_carlo/include/trial_factory.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/profile_cpu.h"

namespace feasst {

FEASST_MAPPER(ProfileCPU,);

ProfileCPU::ProfileCPU(argtype * args) : AnalyzeWriteOnly(args) {}
ProfileCPU::ProfileCPU(argtype args) : ProfileCPU(&args) {
  feasst_check_all_used(args);
}

void ProfileCPU::initialize(MonteCarlo * mc) {
  Analyze::initialize(mc);
  mc->set_timer();
  mc->get_trial_factory()->set_timer();
  mc->get_analyze_factory()->set_timer();
  mc->get_modify_factory()->set_timer();
  printer(header(*mc), output_file(mc->criteria()));
}

std::string ProfileCPU::header(const MonteCarlo& mc) const {
  std::stringstream ss;
  ss << "#percentage of CPU utilized for the following:" << std::endl;
  if (mc.timer()) {
    ss << "all trials,all analyze,all modify,checkpoint,other,,";
  }
  const TimerRDTSC * trial_timer = mc.trials().timer();
  if (trial_timer) {
    const std::vector<double>& trial_times = trial_timer->percents();
    for (const std::shared_ptr<Trial>& trial : mc.trials().trials()) {
      ss << trial->name_or_description() << ",";
    }
    ss << ",";
  }
  const TimerRDTSC * analyze_timer = mc.analyze_factory().timer();
  if (analyze_timer) {
    const std::vector<double>& analyze_times = analyze_timer->percents();
    for (const std::shared_ptr<Analyze>& analyze : mc.analyze_factory().analyzers()) {
      ss << analyze->class_name() << ",";
    }
    ss << ",";
  }
  const TimerRDTSC * modify_timer = mc.modify_factory().timer();
  if (modify_timer) {
    const std::vector<double>& modify_times = modify_timer->percents();
    for (const std::shared_ptr<Modify>& modify : mc.modify_factory().modifiers()) {
      ss << modify->class_name() << ",";
    }
    ss << ",";
  }
  ss << std::endl;
  return ss.str();
}

std::string ProfileCPU::write(const MonteCarlo& mc) {
  std::stringstream ss;
  #ifndef IS_X86
    ss << "ProfileCPU only works on x86 architectures." << std::endl;
    return ss.str();
  #endif
  if (rewrite_header()) {
    ss << header(mc);
  }
  std::vector<double> mc_times;
  if (mc.timer()) {
    mc_times = mc.timer()->percents();
    ss << feasst_str(mc_times) << ",";
  } else {
    mc_times = std::vector<double>(5, 100.);
  }
  const TimerRDTSC * trial_timer = mc.trials().timer();
  if (trial_timer) {
    const std::vector<double>& trial_times = trial_timer->percents();
    for (const double time : trial_times) {
      ss << time*mc_times[0]/100. << ",";
    }
    ss << ",";
  }
  const TimerRDTSC * analyze_timer = mc.analyze_factory().timer();
  if (analyze_timer) {
    const std::vector<double>& analyze_times = analyze_timer->percents();
    for (const double time : analyze_times) {
      ss << time*mc_times[1]/100. << ",";
    }
    ss << ",";
  }
  const TimerRDTSC * modify_timer = mc.modify_factory().timer();
  if (modify_timer) {
    const std::vector<double>& modify_times = modify_timer->percents();
    for (const double time : modify_times) {
      ss << time*mc_times[2]/100. << ",";
    }
  }
  ss << std::endl;
  return ss.str();
}

void ProfileCPU::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(7176, ostr);
}

ProfileCPU::ProfileCPU(std::istream& istr) : AnalyzeWriteOnly(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 7176, "mismatch version:" << version);
}

}  // namespace feasst
