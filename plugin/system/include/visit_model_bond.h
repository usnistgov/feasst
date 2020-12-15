
#ifndef FEASST_SYSTEM_VISIT_MODEL_BOND_H_
#define FEASST_SYSTEM_VISIT_MODEL_BOND_H_

#include <memory>
#include "utils/include/arguments.h"
#include "system/include/visit_model.h"

namespace feasst {

/**
  Interactions between bonded sites using "inter"-like models are computed.
 */
class VisitModelBond : public VisitModel {
 public:
  explicit VisitModelBond(const argtype& args = argtype()) {
    class_name_ = "VisitModelBond"; }
  void compute(
      ModelTwoBody * model,
      const ModelParams& model_params,
      const Select& selection,
      Configuration * config,
      const int group_index) override;
  void compute(
      ModelTwoBody * model,
      const ModelParams& model_params,
      Configuration * config,
      const int group_index) override;
  std::shared_ptr<VisitModel> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit VisitModelBond(std::istream& istr);
  virtual ~VisitModelBond() {}
};

inline std::shared_ptr<VisitModelBond> MakeVisitModelBond(
    const argtype& args = argtype()) {
  return std::make_shared<VisitModelBond>(args);
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_VISIT_MODEL_BOND_H_
