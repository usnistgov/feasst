#include <sstream>
#include "utils/include/serialize.h"
#include "example/include/analyze_example.h"

namespace feasst {

class MapAnalyzeExample {
 public:
  MapAnalyzeExample() {
    auto obj = MakeAnalyzeExample();
    obj->deserialize_map()["AnalyzeExample"] = obj;
  }
};

static MapAnalyzeExample mapper_ = MapAnalyzeExample();

AnalyzeExample::AnalyzeExample(argtype * args) : Analyze(args) {
  group_index_ = integer("group_index", args, 0);
}
AnalyzeExample::AnalyzeExample(argtype args) : AnalyzeExample(&args) {
  FEASST_CHECK_ALL_USED(args);
}

std::string AnalyzeExample::header(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) const {
  std::stringstream ss;
  ss << "dim," << center_[0].status_header() << std::endl;
  return ss.str();
}

void AnalyzeExample::initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  center_.resize(system->configuration().dimension(), feasst::Accumulator());
}

/*
  Create a derived class of LoopConfigOneBody that loops through each site
  and allows some operation for each site by override of the work method.
 */
class AveragePosition : public LoopConfigOneBody {
 public:
  AveragePosition(Position * average) { average_ = average; }
  void work(const Site& site,
    const Configuration& config,
    const LoopDescriptor& data) override { average_->add(site.position()); }
 private:
  Position * average_;
};

void AnalyzeExample::update(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  Position average(system.configuration().dimension());
  AveragePosition loop(&average);
  VisitConfiguration().loop(system.configuration(), &loop, group_index_);
  for (int dim = 0; dim < system.configuration().dimension(); ++dim) {
    center_[dim].accumulate(average.coord(dim));
  }
}

std::string AnalyzeExample::write(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  std::stringstream ss;
  ss << header(criteria, system, trial_factory);
  for (int dim = 0; dim < system.configuration().dimension(); ++dim) {
    ss << dim << "," << center_[dim].status() << std::endl;
  }
  return ss.str();
}

void AnalyzeExample::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(1609, ostr);
  feasst_serialize(group_index_, ostr);
  feasst_serialize_fstobj(center_, ostr);
}

AnalyzeExample::AnalyzeExample(std::istream& istr) : Analyze(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 1609, "version mismatch:" << version);
  feasst_deserialize(&group_index_, istr);
  feasst_deserialize_fstobj(&center_, istr);
}

}  // namespace feasst
