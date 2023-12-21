
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
  Similar to VisitModel, except that the distance between the first sites on
  two different particles are compared against the cutoff + 2*cutoff_outer
  (given in the fstprt files),
  and if that distance is beyond this outer cutoff, then all other interactions
  between these particles.
  Thus, cutoff_outer should be set to the maximum possible distance between
  the first site of the particle and any other site in the same particle.
 */
class VisitModelCutoffOuter : public VisitModel {
 public:
  /**
    args:
    - energy_cutoff: as described in VisitModel.
   */
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

  // compute interactions between particles in the selection
  void compute_between_selection(
    ModelTwoBody * model,
    const ModelParams& model_params,
    const Select& selection,
    Configuration * config,
    const bool is_old_config,
    Position * relative,
    Position * pbc);

  std::shared_ptr<VisitModel> create(std::istream& istr) const override {
    return std::make_shared<VisitModelCutoffOuter>(istr); }
  std::shared_ptr<VisitModel> create(argtype * args) const override {
    return std::make_shared<VisitModelCutoffOuter>(args); }
  void serialize(std::ostream& ostr) const override;
  explicit VisitModelCutoffOuter(std::istream& istr);
  virtual ~VisitModelCutoffOuter() {}

 private:
  double energy_cutoff_;
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
