#include <cmath>
#include "utils/include/utils.h"
#include "utils/include/serialize.h"
#include "math/include/utils_math.h"
#include "configuration/include/domain.h"
#include "system/include/ideal_gas.h"
#include "steppers/include/pair_distribution.h"

namespace feasst {

class MapPairDistributionInner {
 public:
  MapPairDistributionInner() {
    PairDistributionInner().deserialize_map()["PairDistributionInner"] =
      std::make_shared<PairDistributionInner>();
  }
};

static MapPairDistributionInner mapper_inner_ = MapPairDistributionInner();

PairDistributionInner::PairDistributionInner() {
  class_name_ = "PairDistributionInner";
}

void PairDistributionInner::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_pair_distribution_inner_(ostr);
}

void PairDistributionInner::serialize_pair_distribution_inner_(std::ostream& ostr) const {
  feasst_serialize_version(2947, ostr);
  feasst_serialize_fstobj(radial_, ostr);
}

PairDistributionInner::PairDistributionInner(std::istream& istr) : ModelTwoBody(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(2947 == version, version);
  feasst_deserialize_fstobj(&radial_, istr);
}

double PairDistributionInner::energy(
    const double squared_distance,
    const int type1,
    const int type2,
    const ModelParams& model_params) {
  //DEBUG(squared_distance << " " << type1 << " " << type2);
  const double distance = std::sqrt(squared_distance);
  radial_[type1][type2].add(distance);
  radial_[type2][type1].add(distance);
  return 0.;
}

class MapPairDistribution {
 public:
  MapPairDistribution() {
    auto obj = MakePairDistribution();
    obj->deserialize_map()["PairDistribution"] = obj;
  }
};

static MapPairDistribution mapper_ = MapPairDistribution();

PairDistribution::PairDistribution(argtype args) : Modify(&args) {
  //DEBUG("class " << class_name());
  //DEBUG("args " << args_.status());
  dr_ = dble("dr", &args, 0.1);
  DEBUG("is multistate? " << is_multistate());
  DEBUG("is multistate agg? " << is_multistate_aggregate());
  DEBUG("file name " << file_name());
  check_all_used(args);
}

void PairDistribution::initialize(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  DEBUG("init");
  num_updates_ = 0;
  const int num_site_types = system->configuration().num_site_types();
  inter_.radial_.clear();
  intra_.radial_.clear();
  resize(num_site_types, num_site_types, &inter_.radial_);
  resize(num_site_types, num_site_types, &intra_.radial_);
  for (int itype = 0; itype < num_site_types; ++itype) {
    for (int jtype = 0; jtype < num_site_types; ++jtype) {
      Histogram hist;
      hist.set_width_center(dr_, 0.5*dr_);
      inter_.radial_[itype][jtype] = hist;
      intra_.radial_[itype][jtype] = hist;
    }
  }

  // set cutoff to half the minimum box length
  params_ = ModelParams(system->configuration().model_params());
  //DEBUG(params_.cutoff().size());
  const double min_side = system->configuration().domain().min_side_length();
  for (int itype = 0; itype < num_site_types; ++itype) {
    params_.set("cutoff", itype, 0.5*min_side);
  }
  //DEBUG(params_.cutoff().size());
  //DEBUG(params_.cutoff().mixed_value(0, 0));
}

std::string PairDistribution::header(const Criteria& criteria,
    const System& system,
    const TrialFactory& trial_factory) const {
  std::stringstream ss;
  ss << "r,";
  const int num_site_types = system.configuration().num_site_types();
  for (int itype = 0; itype < num_site_types; ++itype) {
    for (int jtype = itype; jtype < num_site_types; ++jtype) {
      ss << itype << "-" << jtype << ",";
    }
  }
  ss << std::endl;
  return ss.str();
}

void PairDistribution::update(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  DEBUG("updating");
  ++num_updates_;
  //DEBUG(params_.cutoff().size());
  //DEBUG(params_.cutoff().mixed_value(0, 0));
  inter_visit_.compute(&inter_, params_, system->get_configuration(), 0);
  intra_visit_.compute(&intra_, params_, system->get_configuration(), 0);
}

std::string PairDistribution::write(Criteria * criteria,
    System * system,
    TrialFactory * trial_factory) {
  std::stringstream ss;
  const grtype& rad = radial(system->configuration());
  for (const grbintype& gr : rad) {
    ss << gr.first << ",";
    for (int itype = 0; itype < static_cast<int>(gr.second.size()); ++itype) {
      for (int jtype = itype; jtype < static_cast<int>(gr.second.size()); ++jtype) {
        ss << gr.second[itype][jtype] << ",";
      }
    }
    ss << std::endl;
  }
  return ss.str();
}

const grtype& PairDistribution::radial(const Configuration& config) {
  ASSERT(num_updates_ > 0, "cannot obtain radial with no updates");
  const int num_site_types = config.num_site_types();
  const std::vector<int> num_sites_of_type = config.num_sites_of_type();
  const int max_bin = static_cast<int>(0.5*config.domain().min_side_length()/dr_);
  for (int bin = 0; bin < max_bin && bin < inter_.radial_[0][0].size(); ++bin) {
    DEBUG("bin: " << bin << " of " << max_bin);
    const double distance = inter_.radial_[0][0].center_of_bin(bin);
    const double lower = inter_.radial_[0][0].edges()[bin];
    const double upper = inter_.radial_[0][0].edges()[bin + 1];
    const double dv = spherical_shell_volume(lower, upper, config.dimension())/
                      config.domain().volume();
    if (static_cast<int>(radial_.size()) <= bin) {
      radial_.push_back(grbintype());
      resize(num_site_types, num_site_types, &radial_[bin].second);
    }
    radial_[bin].first = distance;
    for (int itype = 0; itype < num_site_types; ++itype) {
      const double num_itype = num_sites_of_type[itype];
      for (int jtype = 0; jtype < num_site_types; ++jtype) {
        double norm_fac = 0.;
        if (itype == jtype) norm_fac = 1.;
        const double num_jtype = num_sites_of_type[jtype];
        const Histogram& hist = inter_.radial_[itype][jtype];
        double grbin = 0;
        if (bin < hist.size()) {
          grbin = hist.histogram()[bin]
            /(num_itype - norm_fac)
            /num_jtype
            /static_cast<double>(num_updates_)
            /dv;
        }
        radial_[bin].second[itype][jtype] = grbin;
      }
    }
  }
  return radial_;
}

void PairDistribution::serialize(std::ostream& ostr) const {
  Stepper::serialize(ostr);
  feasst_serialize_version(2034, ostr);
  feasst_serialize(dr_, ostr);
  feasst_serialize_fstobj(inter_visit_, ostr);
  feasst_serialize_fstobj(inter_, ostr);
  feasst_serialize_fstobj(intra_visit_, ostr);
  feasst_serialize_fstobj(intra_, ostr);
  feasst_serialize_fstobj(params_, ostr);
  feasst_serialize(num_updates_, ostr);
}

PairDistribution::PairDistribution(std::istream& istr) : Modify(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(version == 2034, "mismatch version:" << version);
  feasst_deserialize(&dr_, istr);
  feasst_deserialize_fstobj(&inter_visit_, istr);
  feasst_deserialize_fstobj(&inter_, istr);
  feasst_deserialize_fstobj(&intra_visit_, istr);
  feasst_deserialize_fstobj(&intra_, istr);
  feasst_deserialize_fstobj(&params_, istr);
  feasst_deserialize(&num_updates_, istr);
}

PairDistribution::PairDistribution(const Modify& pair_distribution) {
  std::stringstream ss;
  pair_distribution.serialize(ss);
  *this = PairDistribution(ss);
}

}  // namespace feasst
