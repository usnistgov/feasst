
#ifndef FEASST_SYSTEM_LONG_RANGE_CORRECTIONS_H_
#define FEASST_SYSTEM_LONG_RANGE_CORRECTIONS_H_

#include <vector>
#include <sstream>
#include <algorithm>
#include "configuration/include/configuration.h"
#include "system/include/visit_model.h"

namespace feasst {

// HWH: determining number of sites of type is inefficient (order N)
/**
  See Allen and Tildesley or Frenkel and Smit.
 */
class LongRangeCorrections : public VisitModel {
 public:
  LongRangeCorrections() {}

//  // compute number of sites of each type in selection
//  std::vector<int> types(const Select& selection, const Configuration * config) {
//    std::vector<int> count(config->num_site_types());
//    for (int select_index = 0;
//         select_index < selection.num_particles();
//         ++select_index) {
//      const int part_index = selection.particle_index(select_index);
//      const Particle& part = config->select_particle(part_index);
//      for (int site_index : selection.site_indices(select_index)) {
//        const Site& site = part.site(site_index);
//        ++count[site.type()];
//      }
//    }
//    return count;
//  }

  void compute(
      const ModelOneBody& model,
      const ModelParams& model_params,
      const Select& selection,
      Configuration * config,
      const int group_index) override {
    const std::vector<int> num_of_site_type =
      config->num_sites_of_type(group_index);
    std::vector<int> select_types = config->num_sites_of_type(selection);
    DEBUG("sel: " << selection.str());
    DEBUG("num sites of type in selection: " << feasst_str(select_types));
    double en = 0.;
    for (int type1 = 0; type1 < config->num_site_types(); ++type1) {
      const double num_type1 = num_of_site_type[type1];
      const double num_type1_sel = select_types[type1];
      for (int type2 = 0; type2 < config->num_site_types(); ++type2) {
        const double num_type2 = num_of_site_type[type2];
        const double num_type2_sel = select_types[type2];
        en += (num_type1*num_type2_sel + num_type1_sel*num_type2 - num_type1_sel*num_type2_sel)
          *energy_(type1, type2, config, model_params);
      }
    }
    set_energy(en);
  }

  void compute(
      const ModelOneBody& model,
      const ModelParams& model_params,
      Configuration * config,
      const int group_index = 0) override {
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

  void serialize(std::ostream& ostr) const override {
    ostr << class_name_ << " ";
    serialize_visit_model_(ostr);
    feasst_serialize_version(874, ostr);
  }

  std::shared_ptr<VisitModel> create(std::istream& istr) const override {
    return std::make_shared<LongRangeCorrections>(istr);
  }

  LongRangeCorrections(std::istream& istr) : VisitModel(istr) {
    const int version = feasst_deserialize_version(istr);
    ASSERT(874 == version, version);
  }

 private:
  const std::string class_name_ = "LongRangeCorrections";

  double energy_(
      const int type1,
      const int type2,
      const Configuration * config,
      const ModelParams& model_params) const {
    const double epsilon = model_params.mixed_epsilon()[type1][type2];
    const double sigma = model_params.mixed_sigma()[type1][type2];
    const double cutoff = model_params.mixed_cutoff()[type1][type2];
    const double prefactor = epsilon*(8./3.)*PI*pow(sigma, 3)*
      ((1./3.)*pow(sigma/cutoff, 9) - pow(sigma/cutoff, 3));
    return prefactor/config->domain()->volume();
  }
};

inline std::shared_ptr<LongRangeCorrections> MakeLongRangeCorrections() {
  return std::make_shared<LongRangeCorrections>();
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_LONG_RANGE_CORRECTIONS_H_
