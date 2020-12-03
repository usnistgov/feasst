
#ifndef OPT_RPM_VISIT_MODEL_OPT_RPM_H_
#define OPT_RPM_VISIT_MODEL_OPT_RPM_H_

#include "system/include/visit_model.h"

namespace feasst {

/**
  Note: Do not use this, it is still being implemented.
  Optimize the visitor for single-site RPM interactions in 3D.
 */
class VisitModelOptRPM : public VisitModel {
 public:
  VisitModelOptRPM() : VisitModel() {
    class_name_ = "VisitModelOptRPM"; }
  VisitModelOptRPM(std::shared_ptr<VisitModelInner> inner) : VisitModel(inner) {}
  void precompute(Configuration * config) override;
  void compute(
      ModelTwoBody * model,
      const ModelParams& model_params,
      const Select& selection,
      Configuration * config,
      const int group_index = 0) override;

  std::shared_ptr<VisitModel> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  VisitModelOptRPM(std::istream& istr);
  ~VisitModelOptRPM() {}

 private:
  double alpha_;
  double conversion_factor_;
};

inline std::shared_ptr<VisitModelOptRPM> MakeVisitModelOptRPM() {
  return std::make_shared<VisitModelOptRPM>();
}

inline std::shared_ptr<VisitModelOptRPM> MakeVisitModelOptRPM(
    std::shared_ptr<VisitModelInner> inner) {
  return std::make_shared<VisitModelOptRPM>(inner);
}

}  // namespace feasst

#endif  // OPT_RPM_VISIT_MODEL_OPT_RPM_H_
