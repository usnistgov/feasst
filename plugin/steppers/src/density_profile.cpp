#include <deque>
#include <cmath>
#include "utils/include/utils.h"
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "configuration/include/domain.h"
#include "configuration/include/visit_configuration.h"
#include "system/include/ideal_gas.h"
#include "steppers/include/density_profile.h"

namespace feasst {

class MapDensityProfile {
 public:
  MapDensityProfile() {
    auto obj = MakeDensityProfile();
    obj->deserialize_map()["DensityProfile"] = obj;
  }
};

static MapDensityProfile mapper_ = MapDensityProfile();

DensityProfile::DensityProfile(argtype args) : Analyze(&args) {
  dimension_ = integer("dimension", &args, 0);
  dr_ = dble("dr", &args, 0.1);
  center_ = dble("center", &args, 0.);
  check_all_used(args);
}

void DensityProfile::initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  DEBUG("init");
  const int num_site_types = system->configuration().num_site_types();
  data_.resize(num_site_types);
  const double max_side = system->configuration().domain().max_side_length();
  for (int type = 0; type < num_site_types; ++type) {
    Histogram hist;
    hist.set_width_center(dr_, center_);
    hist.add(-0.5*max_side, false);
    hist.add(0.5*max_side, false);
    data_[type] = hist;
  }
}

std::string DensityProfile::header(const Criteria& criteria,
    const System& system,
    const TrialFactory& trials) const {
  std::stringstream ss;
  ss << "r,";
  for (int type = 0; type < system.configuration().num_site_types(); ++type) {
    ss << type << ",";
  }
  ss << std::endl;
  return ss.str();
}

class ComputeProfile : public LoopConfigOneBody {
 public:
  ComputeProfile(const int dim, std::vector<Histogram> * hist)
    : dim_(dim), hist_(hist) {}
  void work(const Site& site,
      const Configuration& config,
      const LoopDescriptor& data) override {
    (*hist_)[site.type()].add(site.position().coord(dim_));
  }
 private:
  int dim_;
  std::vector<Histogram> * hist_;
};

void DensityProfile::update(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  DEBUG("updating");
  ComputeProfile comp(dimension_, &data_);
  VisitConfiguration().loop(system.configuration(), &comp);
}

std::vector<std::vector<std::vector<double> > > DensityProfile::profile() const {
  std::vector<std::vector<std::vector<double> > > prof;
  resize(data_[0].size(), data_.size(), 2, &prof);
  std::vector<double> sum;
  for (int type = 0; type < static_cast<int>(data_.size()); ++type) {
    const std::deque<double>& data = data_[type].histogram();
    sum.push_back(std::accumulate(data.begin(), data.end(), 0.));
  }
  for (int bin = 0; bin < data_[0].size(); ++bin) {
    for (int type = 0; type < static_cast<int>(data_.size()); ++type) {
      prof[bin][type][0] = data_[type].center_of_bin(bin);
      prof[bin][type][1] = data_[type].histogram()[bin]/sum[type];
    }
  }
  return prof;
}

std::string DensityProfile::write(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) {
  std::stringstream ss;
  ss << header(criteria, system, trial_factory);
  const std::vector<std::vector<std::vector<double> > > prof = profile();
  for (const std::vector<std::vector<double> >& profbin : prof) {
    for (int type = 0; type < static_cast<int>(profbin.size()); ++type) {
      const std::vector<double>& profbintype = profbin[type];
      if (type == 0) ss << profbintype[0] << ",";
      ss << profbintype[1] << ",";
    }
    ss << std::endl;
  }
  return ss.str();
}

void DensityProfile::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(9687, ostr);
  feasst_serialize(dimension_, ostr);
  feasst_serialize(dr_, ostr);
  feasst_serialize(center_, ostr);
  feasst_serialize_fstobj(data_, ostr);
}

DensityProfile::DensityProfile(std::istream& istr) : Analyze(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 9687, "mismatch version:" << version);
  feasst_deserialize(&dimension_, istr);
  feasst_deserialize(&dr_, istr);
  feasst_deserialize(&center_, istr);
  feasst_deserialize_fstobj(&data_, istr);
}

DensityProfile::DensityProfile(const Analyze& density_profile) {
  std::stringstream ss;
  density_profile.serialize(ss);
  *this = DensityProfile(ss);
}

}  // namespace feasst
