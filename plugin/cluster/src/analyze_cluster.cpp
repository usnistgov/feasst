#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "math/include/random_mt19937.h"
#include "monte_carlo/include/trial_select.h"
#include "monte_carlo/include/trial_factory.h"
#include "monte_carlo/include/trial_stage.h"
#include "cluster/include/analyze_cluster.h"

namespace feasst {

FEASST_MAPPER(AnalyzeCluster,);

AnalyzeCluster::AnalyzeCluster(argtype * args) : Analyze(args) {}
AnalyzeCluster::AnalyzeCluster(argtype args) : AnalyzeCluster(&args) {
  feasst_check_all_used(args);
}

void AnalyzeCluster::initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  printer(header(*criteria, *system, *trial_factory),
          output_file(*criteria));
}

std::string AnalyzeCluster::header(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) const {
  std::stringstream ss;
  ss << accumulator().status_header() << std::endl;
  return ss.str();
}

void AnalyzeCluster::update(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  const int trial = trial_factory.last_index();
  ASSERT(trial != -1, "no trials to harvest cluster information from.");
  const TrialSelect& sel = trial_factory.trial(trial).stage(0).trial_select();
  if (sel.class_name() == "SelectCluster") {
    const Accumulator& csize = sel.printable("cluster_size");
    if (csize.num_values() > 0) {
      get_accumulator()->accumulate(csize.last_value());
    }
  }
}

std::string AnalyzeCluster::write(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  std::stringstream ss;
  if (rewrite_header()) {
    ss << header(criteria, system, trial_factory);
  }
  ss << accumulator().status() << std::endl;
  DEBUG(ss.str());
  return ss.str();
}

void AnalyzeCluster::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(325, ostr);
}

AnalyzeCluster::AnalyzeCluster(std::istream& istr) : Analyze(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 325, "mismatch version:" << version);
}

AnalyzeCluster::AnalyzeCluster(const Analyze& energy) {
  std::stringstream ss;
  energy.serialize(ss);
  *this = AnalyzeCluster(ss);
}

}  // namespace feasst
