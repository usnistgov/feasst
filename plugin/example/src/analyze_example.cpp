#include "utils/include/arguments.h"
#include "utils/include/serialize.h"
#include "math/include/accumulator.h"
#include "configuration/include/particle_factory.h"
#include "configuration/include/configuration.h"
#include "configuration/include/visit_configuration.h"
#include "system/include/system.h"
#include "monte_carlo/include/monte_carlo.h"
#include "example/include/analyze_example.h"

namespace feasst {

FEASST_MAPPER(AnalyzeExample,);

AnalyzeExample::AnalyzeExample(argtype * args) : Analyze(args) {
  group_ = str("group", args, "");
}

AnalyzeExample::AnalyzeExample(argtype args) : AnalyzeExample(&args) {
  feasst_check_all_used(args);
}

AnalyzeExample::~AnalyzeExample() {}

std::string AnalyzeExample::header(const MonteCarlo& mc) const {
  std::stringstream ss;
  ss << "dim," << center_[0]->status_header() << std::endl;
  return ss.str();
}

void AnalyzeExample::initialize(MonteCarlo * mc) {
  Analyze::initialize(mc);
  if (!group_.empty()) {
    group_index_ = mc->configuration().group_index(group_);
  }
  center_.resize(configuration(mc->system()).dimension());
  for (std::unique_ptr<Accumulator>& acc : center_) {
    acc = std::make_unique<Accumulator>();
  }
}

/*
  Create a derived class of LoopConfigOneBody that loops through each site
  and allows some operation for each site by override of the work method.
 */
class AveragePosition : public LoopConfigOneBody {
 public:
  AveragePosition() {}
  void work(const Site& site,
    const Configuration& config,
    const LoopDescriptor& data) override { average_.add(site.position()); }
  const Position& average() const { return average_; }
  void zero_average(const int num_dimen) { average_.set_to_origin(num_dimen); }
 private:
  Position average_;
};

void AnalyzeExample::update(const MonteCarlo& mc) {
  const Configuration& config = configuration(mc.system());
  if (!loop_) {
    loop_ = std::make_unique<AveragePosition>();
    visit_ = std::make_unique<VisitConfiguration>();
  }
  loop_->zero_average(config.dimension());
  visit_->loop(config, loop_.get(), group_index_);
  if (visit_->num_sites() > 0) {
    for (int dim = 0; dim < config.dimension(); ++dim) {
      center_[dim]->accumulate(loop_->average().coord(dim));
    }
  }
}

std::string AnalyzeExample::write(const MonteCarlo& mc) {
  std::stringstream ss;
  ss << header(mc);
  for (int dim = 0; dim < configuration(mc.system()).dimension(); ++dim) {
    ss << dim << "," << center_[dim]->status() << std::endl;
  }
  return ss.str();
}

const Accumulator& AnalyzeExample::geometric_center(const int dimension) const {
  return *center_[dimension];
}

const std::vector<std::unique_ptr<Accumulator> >&
  AnalyzeExample::geometric_center() const { return center_; }

void AnalyzeExample::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(1609, ostr);
  feasst_serialize(group_index_, ostr);
  feasst_serialize(center_, ostr);
}

AnalyzeExample::AnalyzeExample(std::istream& istr) : Analyze(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 1609, "version mismatch:" << version);
  feasst_deserialize(&group_index_, istr);
  feasst_deserialize(&center_, istr);
}

}  // namespace feasst
