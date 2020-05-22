#include <cmath>
#include "system/include/long_range_corrections.h"
#include "utils/include/serialize.h"
#include "math/include/constants.h"

namespace feasst {

class MapLongRangeCorrections {
 public:
  MapLongRangeCorrections() {
    LongRangeCorrections().deserialize_map()["LongRangeCorrections"] =
      MakeLongRangeCorrections();
  }
};

static MapLongRangeCorrections mapper_ = MapLongRangeCorrections();

void LongRangeCorrections::serialize(std::ostream& ostr) const {
  ostr << class_name_ << " ";
  serialize_visit_model_(ostr);
  feasst_serialize_version(874, ostr);
}

LongRangeCorrections::LongRangeCorrections(std::istream& istr)
  : VisitModel(istr) {
  const int version = feasst_deserialize_version(istr);
  ASSERT(874 == version, version);
}

void LongRangeCorrections::compute(
    const ModelOneBody& model,
    const ModelParams& model_params,
    const Select& selection,
    Configuration * config,
    const int group_index) {
  const std::vector<int> num_of_site_type =
    config->num_sites_of_type(group_index);
  std::vector<int> select_types = config->num_sites_of_type(selection);
  DEBUG("sel: " << selection.str());
  DEBUG("num sites of type in selection: " << feasst_str(select_types));
  double en = 0.;
  for (int type1 = 0; type1 < config->num_site_types(); ++type1) {
    DEBUG("type1 " << type1);
    const double num_type1 = num_of_site_type[type1];
    const double num_type1_sel = select_types[type1];
    DEBUG("num_particles " << config->num_particles());
    DEBUG("num_type1 " << num_type1);
    DEBUG("num_type1_sel " << num_type1_sel);
    for (int type2 = 0; type2 < config->num_site_types(); ++type2) {
      DEBUG("type2 " << type2);
      const double num_type2 = num_of_site_type[type2];
      const double num_type2_sel = select_types[type2];
      DEBUG("num_type2 " << num_type2);
      DEBUG("num_type2_sel " << num_type2_sel);
      en += (num_type1*num_type2_sel +
             num_type1_sel*num_type2 -
             num_type1_sel*num_type2_sel)
        *energy_(type1, type2, config, model_params);
    }
  }
  set_energy(en);
}

void LongRangeCorrections::compute(
    const ModelOneBody& model,
    const ModelParams& model_params,
    Configuration * config,
    const int group_index) {
  double en = 0;
  const std::vector<int> num_of_site_type =
    config->num_sites_of_type(group_index);
  DEBUG("num sites of type in group: " << feasst_str(num_of_site_type));
  for (int type1 = 0; type1 < config->num_site_types(); ++type1) {
    for (int type2 = 0; type2 < config->num_site_types(); ++type2) {
      en += num_of_site_type[type1]*num_of_site_type[type2]*
        energy_(type1, type2, config, model_params);
    }
  }
  set_energy(en);
}

double LongRangeCorrections::energy_(
    const int type1,
    const int type2,
    const Configuration * config,
    const ModelParams& model_params) const {
  const double epsilon = model_params.mixed_epsilon()[type1][type2];
  const double sigma = model_params.mixed_sigma()[type1][type2];
  const double cutoff = model_params.mixed_cutoff()[type1][type2];
  const double prefactor = epsilon*(8./3.)*PI*std::pow(sigma, 3)*
    ((1./3.)*std::pow(sigma/cutoff, 9) - std::pow(sigma/cutoff, 3));
  return prefactor/config->domain().volume();
}

}  // namespace feasst
