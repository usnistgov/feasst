#include "utils/include/serialize.h"
#include "utils/include/arguments.h"
#include "math/include/random_mt19937.h"
#include "monte_carlo/include/trial_select.h"
#include "monte_carlo/include/trial_factory.h"
#include "monte_carlo/include/trial_stage.h"
#include "monte_carlo/include/monte_carlo.h"
#include "cluster/include/analyze_cluster.h"

namespace feasst {

FEASST_MAPPER(AnalyzeCluster,);

AnalyzeCluster::AnalyzeCluster(argtype * args) : Analyze(args) {}
AnalyzeCluster::AnalyzeCluster(argtype args) : AnalyzeCluster(&args) {
  feasst_check_all_used(args);
}

void AnalyzeCluster::initialize(MonteCarlo * mc) {
  printer(header(*mc), output_file(mc->criteria()));
}

std::string AnalyzeCluster::header(const MonteCarlo& mc) const {
  std::stringstream ss;
  ss << accumulator().status_header() << std::endl;
  return ss.str();
}

void AnalyzeCluster::update(const MonteCarlo& mc) {
  const TrialFactory& trial_factory = mc.trial_factory();
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

std::string AnalyzeCluster::write(const MonteCarlo& mc) {
  std::stringstream ss;
  if (rewrite_header()) {
    ss << header(mc);
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
