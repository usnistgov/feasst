
#ifndef FEASST_SYSTEM_LONG_RANGE_CORRECTIONS_H_
#define FEASST_SYSTEM_LONG_RANGE_CORRECTIONS_H_

#include <memory>
#include <vector>
#include <string>
#include <sstream>
#include "system/include/visit_model.h"

namespace feasst {

class Configuration;

// HWH: determining number of sites of type is inefficient (order N)
/**
  These are the long range corrections assuming a 12-6 Lennard-Jones potential
  and that the radial distribution function is one beyond the cutoff.

  See Allen and Tildesley or Frenkel and Smit.

  \f$ U_{LRC} = \sum_i \sum_j U_{LRC}^{ij} = \sum_i \sum_j C_{ij} n_i n_j\f$

  where \f$i\f$ and \f$j\f$ are particle types i and j, \f$n\f$ are the number
  of sites of the given type, and the constant, C_{ij} is given by

  \f$C_{ij} = \frac{8\epsilon_{ij}\pi\sigma_{ij}^3}{3 V}\left[\frac{1}{3}\left(\frac{\sigma_{ij}}{r^c_{ij}}\right)^9 - \left(\frac{\sigma_{ij}}{r^c_{ij}}\right)^3\right]\f$

  which depends only on the LennardJones parameters.

  Only trials that change the number of sites of a type contribute to changes
  in energy.

  If a selection of sites are to be deleted, the energetic contribution of the
  selection may be computed as:

  \f$U_{with} - U_{without} \propto n_i n_j - (n_i-n_i^s)(n_j-n_j^2)\f$

  \f$U_{with} - U_{without} \propto n_i^s n_j + n_i n_j^s - n_i^s n_j^s\f$

  where the superscript s refers to the number of sites in the selection.

  If a selection is to be added, the energetic contribution of the selection
  may be computed as:

  \f$U_{with} - U_{without} \propto (n_i+n_i^s)(n_j+n_j^s) - n_i n_j\f$

  \f$U_{with} - U_{without} \propto n_i^s n_j + n_i n_j^s + n_i^s n_j^s\f$
 */
class LongRangeCorrections : public VisitModel {
 public:
  LongRangeCorrections() { class_name_ = "LongRangeCorrections"; }
  void precompute(Configuration * config) override;
  void compute(
      ModelOneBody * model,
      const ModelParams& model_params,
      const Select& selection,
      Configuration * config,
      const int group_index) override;
  void compute(
      ModelOneBody * model,
      const ModelParams& model_params,
      Configuration * config,
      const int group_index = 0) override;
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<VisitModel> create(std::istream& istr) const override {
    return std::make_shared<LongRangeCorrections>(istr);
  }
  std::shared_ptr<VisitModel> create(argtype * args) const override {
    return std::make_shared<LongRangeCorrections>();
  }
  explicit LongRangeCorrections(std::istream& istr);

 private:
  // temporary, and not serialized
  std::vector<int> num_of_site_type_;
  std::vector<int> select_types_;

  double energy_(
    const int type1,
    const int type2,
    const Configuration * config,
    const ModelParams& model_params) const;
};

inline std::shared_ptr<LongRangeCorrections> MakeLongRangeCorrections() {
  return std::make_shared<LongRangeCorrections>();
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_LONG_RANGE_CORRECTIONS_H_
