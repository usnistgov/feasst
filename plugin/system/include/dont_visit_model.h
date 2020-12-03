
#ifndef FEASST_SYSTEM_DONT_VISIT_MODEL_H_
#define FEASST_SYSTEM_DONT_VISIT_MODEL_H_

#include <memory>
#include "system/include/visit_model.h"

namespace feasst {

/**
  Return zero energy instead of visiting a model.
  This may be used for reference potentials of single step trials.
 */
class DontVisitModel : public VisitModel {
 public:
  DontVisitModel() : VisitModel() {
    class_name_ = "DontVisitModel"; }
  void compute(
      ModelOneBody * model,
      const ModelParams& model_params,
      Configuration * config,
      const int group_index = 0) override { set_energy(0.); }
  void compute(
      ModelOneBody * model,
      const ModelParams& model_params,
      const Select& selection,
      Configuration * config,
      const int group_index = 0) override { set_energy(0.); }
  void compute(
      ModelTwoBody * model,
      const ModelParams& model_params,
      Configuration * config,
      const int group_index = 0) override { set_energy(0.); }
  void compute(
      ModelTwoBody * model,
      const ModelParams& model_params,
      const Select& selection,
      Configuration * config,
      const int group_index = 0) override { set_energy(0.); }
  void compute(
      ModelThreeBody * model,
      const ModelParams& model_params,
      Configuration * config,
      const int group_index = 0) override { set_energy(0.); }
  std::shared_ptr<VisitModel> create(std::istream& istr) const override;
  void serialize(std::ostream& ostr) const override;
  explicit DontVisitModel(std::istream& istr);
  virtual ~DontVisitModel() {}
};

inline std::shared_ptr<DontVisitModel> MakeDontVisitModel() {
  return std::make_shared<DontVisitModel>();
}

}  // namespace feasst

#endif  // FEASST_SYSTEM_DONT_VISIT_MODEL_H_
