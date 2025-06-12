#include <deque>
#include <numeric>  // accumulate
#include "utils/include/arguments.h"
#include "utils/include/utils.h"
#include "utils/include/serialize.h"
#include "configuration/include/site.h"
#include "configuration/include/domain.h"
#include "configuration/include/visit_configuration.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/density_profile.h"

namespace feasst {

FEASST_MAPPER(DensityProfile,);

DensityProfile::DensityProfile(argtype * args) : Analyze(args) {
  dimension_ = integer("dimension", args, 0);
  dr_ = dble("dr", args, 0.1);
  center_ = dble("center", args, 0.);
}
DensityProfile::DensityProfile(argtype args) : DensityProfile(&args) {
  feasst_check_all_used(args);
}

void DensityProfile::initialize(MonteCarlo * mc) {
  Analyze::initialize(mc);
  const Configuration& config = configuration(mc->system());
  DEBUG("init");
  const int num_site_types = config.num_site_types();
  data_.resize(num_site_types);
  const double max_side = config.domain().max_side_length();
  for (int type = 0; type < num_site_types; ++type) {
    Histogram hist;
    hist.set_width_center(dr_, center_);
    hist.add(-0.5*max_side, false);
    hist.add(0.5*max_side, false);
    data_[type] = hist;
  }
}

std::string DensityProfile::header(const MonteCarlo& mc) const {
  const Configuration& config = configuration(mc.system());
  std::stringstream ss;
  ss << "r,";
  for (int type = 0; type < config.num_site_types(); ++type) {
    ss << type << ",";
  }
  ss << std::endl;
  return ss.str();
}

class ComputeProfile : public LoopConfigOneBody {
 public:
  ComputeProfile(const int dim, std::vector<Histogram> * hist, const Domain * domain)
    : dim_(dim), hist_(hist), domain_(domain) {}
  void work(const Site& site,
      const Configuration& config,
      const LoopDescriptor& data) override {
    domain_->wrap_opt(site.position(), &wrapped_, &scaled_);
    (*hist_)[site.type()].add(wrapped_.coord(dim_));
    //(*hist_)[site.type()].add(site.position().coord(dim_));
  }
 private:
  int dim_;
  std::vector<Histogram> * hist_;
  const Domain * domain_;
  Position wrapped_, scaled_;
};

void DensityProfile::update(const MonteCarlo& mc) {
  DEBUG("updating");
  const Configuration& config = configuration(mc.system());
  ComputeProfile comp(dimension_, &data_, &config.domain());
  VisitConfiguration().loop(config, &comp);
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

std::string DensityProfile::write(const MonteCarlo& mc) {
  std::stringstream ss;
  ss << header(mc);
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
