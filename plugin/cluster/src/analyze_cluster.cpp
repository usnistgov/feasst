#include "utils/include/serialize.h"
#include "cluster/include/analyze_cluster.h"

namespace feasst {

class MapAnalyzeCluster {
 public:
  MapAnalyzeCluster() {
    auto obj = MakeAnalyzeCluster();
    obj->deserialize_map()["AnalyzeCluster"] = obj;
  }
};

static MapAnalyzeCluster mapper_ = MapAnalyzeCluster();

AnalyzeCluster::AnalyzeCluster(argtype * args) : Analyze(args) {}
AnalyzeCluster::AnalyzeCluster(argtype args) : AnalyzeCluster(&args) {
  FEASST_CHECK_ALL_USED(args);
}

void AnalyzeCluster::initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  printer(header(*criteria, *system, *trial_factory),
          file_name(*criteria));
}

std::string AnalyzeCluster::header(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) const {
  std::stringstream ss;
  ss << accumulator_.status_header() << std::endl;
  return ss.str();
}

void AnalyzeCluster::update(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  const int trial = trial_factory.last_index();
  const TrialSelect& sel = trial_factory.trial(trial).stage(0).trial_select();
  if (sel.class_name() == "SelectCluster") {
    accumulator_.accumulate(sel.printable("cluster_size").last_value());
  }
}

std::string AnalyzeCluster::write(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  std::stringstream ss;
  if (rewrite_header()) {
    ss << header(criteria, system, trial_factory);
  }
  ss << accumulator_.status() << std::endl;
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
