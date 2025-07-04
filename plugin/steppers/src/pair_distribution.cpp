#include "utils/include/arguments.h"
#include "utils/include/utils.h"
#include "utils/include/serialize_extra.h"
#include "math/include/utils_math.h"
#include "configuration/include/domain.h"
#include "configuration/include/configuration.h"
#include "system/include/system.h"
#include "system/include/visit_model.h"
#include "monte_carlo/include/monte_carlo.h"
#include "steppers/include/pair_distribution.h"

namespace feasst {

FEASST_MAPPER(PairDistributionInner,);

PairDistributionInner::PairDistributionInner() {
  class_name_ = "PairDistributionInner";
}
PairDistributionInner::~PairDistributionInner() {}

void PairDistributionInner::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_pair_distribution_inner_(ostr);
}

void PairDistributionInner::serialize_pair_distribution_inner_(std::ostream& ostr) const {
  serialize_model_(ostr);
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

FEASST_MAPPER(PairDistribution,);

PairDistribution::PairDistribution(argtype * args) : Modify(args) {
  //DEBUG("class " << class_name());
  //DEBUG("args " << args_.status());
  dr_ = dble("dr", args, 0.1);
  print_intra_ = boolean("print_intra", args, false);
  DEBUG("is multistate? " << is_multistate());
  DEBUG("is multistate agg? " << is_multistate_aggregate());
  DEBUG("output_file " << output_file());
}
PairDistribution::PairDistribution(argtype args) : PairDistribution(&args) {
  feasst_check_all_used(args);
}

void PairDistribution::initialize(MonteCarlo * mc) {
  Modify::initialize(mc);
  Configuration * config = mc->get_system()->get_configuration(configuration_index());
  DEBUG("init");
  num_updates_ = 0;
  const int num_site_types = config->num_site_types();
  DEBUG("num_site_types " << num_site_types);
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
  params_ = deep_copy(config->model_params());
  const double min_side = config->domain().inscribed_sphere_diameter();
  DEBUG("min_side " << min_side);
  for (int itype = 0; itype < num_site_types; ++itype) {
    params_.set("cutoff", itype, 0.5*min_side);
    for (int jtype = 0; jtype < num_site_types; ++jtype) {
      params_.set("cutoff", itype, jtype, 0.5*min_side);
    }
  }
  inter_visit_.precompute(config);
  intra_visit_.precompute(config);
}

std::string PairDistribution::header(const MonteCarlo& mc) const {
  const int num_site_types = configuration(mc.system()).num_site_types();
  std::stringstream ss;
  ss << "#\"num_site_types\":" << num_site_types << "," << std::endl;
  ss << "r,";
  for (int itype = 0; itype < num_site_types; ++itype) {
    for (int jtype = itype; jtype < num_site_types; ++jtype) {
      std::stringstream tt;
      tt << itype << "-" << jtype;
      ss << "g" << tt.str() << "," ;
      if (print_intra_) {
        ss << "h"  << tt.str() << "," ;
      }
    }
  }
  ss << std::endl;
  return ss.str();
}

void PairDistribution::update(MonteCarlo * mc) {
  Configuration * config = mc->get_system()->get_configuration(configuration_index());
  DEBUG("updating");
  ++num_updates_;
  //DEBUG(params_.cutoff().size());
  //DEBUG(params_.cutoff().mixed_value(0, 0));
  inter_visit_.compute(&inter_, params_, config, 0);
  intra_visit_.compute(&intra_, params_, config, 0);
}

std::string PairDistribution::write(MonteCarlo * mc) {
  std::stringstream ss;
  ss << header(*mc);
  const Configuration& config = configuration(mc->system());
  if (num_updates_ <= 0) {
    return ss.str();
  }
  const grtype& rad = radial(config);
  const std::vector<std::vector<Histogram> >& hr = intra_.radial_;
  //ASSERT(static_cast<int>(rad.size()) == hr[0][0].size(), "err");
  //for (const grbintype& gr : rad) {
  for (int bin = 0; bin < static_cast<int>(rad.size()); ++bin) {
    const grbintype& gr = rad[bin];
    ss << gr.first << ",";
    for (int itype = 0; itype < static_cast<int>(gr.second.size()); ++itype) {
      for (int jtype = itype; jtype < static_cast<int>(gr.second.size()); ++jtype) {
        ss << gr.second[itype][jtype] << ",";
        if (print_intra_) {
          int h = 0;
          const Histogram& hrij = hr[itype][jtype];
          if (bin < hrij.size()) {
            h = hrij.histogram()[bin]
              /static_cast<double>(num_updates_)
              /static_cast<double>(config.num_particles());
          }
          ss << h << ",";
        }
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
  std::vector<std::vector<int> > num_sites_of_type_in_particle =
    config.num_site_types_per_particle_type();
  const int max_bin = static_cast<int>(0.5*config.domain().inscribed_sphere_diameter()/dr_);
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
      const int ipart_type = config.site_type_to_particle_type(itype);
      for (int jtype = 0; jtype < num_site_types; ++jtype) {
        const int jpart_type = config.site_type_to_particle_type(jtype);
        double norm_fac = 0.;
        if (ipart_type == jpart_type) {
          norm_fac = num_sites_of_type_in_particle[ipart_type][itype];
        }
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
  feasst_serialize(print_intra_, ostr);
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
  feasst_deserialize(&print_intra_, istr);
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
