
#ifndef FEASST_SYSTEM_LONG_RANGE_CORRECTIONS_H_
#define FEASST_SYSTEM_LONG_RANGE_CORRECTIONS_H_

#include <memory>
#include <vector>
#include <string>
#include <sstream>
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
  void compute(
      const ModelOneBody& model,
      const ModelParams& model_params,
      const Select& selection,
      Configuration * config,
      const int group_index) override;
  void compute(
      const ModelOneBody& model,
      const ModelParams& model_params,
      Configuration * config,
      const int group_index = 0) override;
  void serialize(std::ostream& ostr) const override;
  std::shared_ptr<VisitModel> create(std::istream& istr) const override {
    return std::make_shared<LongRangeCorrections>(istr);
  }
  explicit LongRangeCorrections(std::istream& istr);

 private:
  const std::string class_name_ = "LongRangeCorrections";

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
