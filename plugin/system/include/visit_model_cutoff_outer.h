
#ifndef FEASST_SYSTEM_VISIT_MODEL_CUTOFF_OUTER_H_
#define FEASST_SYSTEM_VISIT_MODEL_CUTOFF_OUTER_H_

#include <memory>
#include "utils/include/arguments.h"
#include "system/include/visit_model.h"
#include "system/include/cells.h"

namespace feasst {

class CutoffOuter : public ModelParam {
 public:
  CutoffOuter() : ModelParam() { class_name_ = "cutoff_outer"; }
  std::shared_ptr<ModelParam> create(std::istream& istr) const override {
    return std::make_shared<CutoffOuter>(istr); }
  void serialize(std::ostream& ostr) const override;
  explicit CutoffOuter(std::istream& istr);
  virtual ~CutoffOuter() {}
};

/**
  Compute many-body inter-particle interactions using a cell list.
 */
class VisitModelCutoffOuter : public VisitModel {
 public:
  explicit VisitModelCutoffOuter(argtype args);
  explicit VisitModelCutoffOuter(argtype * args);

  /// Same as above, but with an inner.
  explicit VisitModelCutoffOuter(std::shared_ptr<VisitModelInner> inner,
    argtype args);

  void compute(
      ModelTwoBody * model,
      const ModelParams& model_params,
      const Select& selection,
      Configuration * config,
      const int group_index) override;

  std::shared_ptr<VisitModel> create(std::istream& istr) const {
    return std::make_shared<VisitModelCutoffOuter>(istr); }
  std::shared_ptr<VisitModel> create(argtype * args) const {
    return std::make_shared<VisitModelCutoffOuter>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit VisitModelCutoffOuter(std::istream& istr);
  virtual ~VisitModelCutoffOuter() {}
};

inline std::shared_ptr<VisitModelCutoffOuter> MakeVisitModelCutoffOuter(
    argtype args = argtype()) {
  return std::make_shared<VisitModelCutoffOuter>(args);
}

inline std::shared_ptr<VisitModelCutoffOuter> MakeVisitModelCutoffOuter(
    std::shared_ptr<VisitModelInner> inner,
    argtype args = argtype()) {
  return std::make_shared<VisitModelCutoffOuter>(inner, args);
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_VISIT_MODEL_CUTOFF_OUTER_H_
