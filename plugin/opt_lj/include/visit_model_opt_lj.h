
#ifndef OPT_LJ_VISIT_MODEL_OPT_LJ_H_
#define OPT_LJ_VISIT_MODEL_OPT_LJ_H_

#include "system/include/visit_model.h"

namespace feasst {

/**
  Optimize the visitor for single-site LJ interactions in 3D.
 */
class VisitModelOptLJ : public VisitModel {
 public:
  VisitModelOptLJ() : VisitModel() {
    class_name_ = "VisitModelOptLJ"; }
  VisitModelOptLJ(std::shared_ptr<VisitModelInner> inner) : VisitModel(inner) {}
  void compute(
      const ModelTwoBody& model,
      const ModelParams& model_params,
      const Select& selection,
      Configuration * config,
      const int group_index = 0) override;

  std::shared_ptr<VisitModel> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  VisitModelOptLJ(std::istream& istr);
  ~VisitModelOptLJ() {}

 private:
  // temporary
  double squared_distance_;
};

inline std::shared_ptr<VisitModelOptLJ> MakeVisitModelOptLJ() {
  return std::make_shared<VisitModelOptLJ>();
}

inline std::shared_ptr<VisitModelOptLJ> MakeVisitModelOptLJ(
    std::shared_ptr<VisitModelInner> inner) {
  return std::make_shared<VisitModelOptLJ>(inner);
}

}  // namespace feasst

#endif  // OPT_LJ_VISIT_MODEL_OPT_LJ_H_
