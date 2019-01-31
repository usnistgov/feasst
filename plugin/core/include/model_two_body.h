
#ifndef FEASST_CORE_MODEL_TWO_BODY_H_
#define FEASST_CORE_MODEL_TWO_BODY_H_

#include "core/include/model.h"
#include "core/include/visit_model.h"
#include "core/include/particle.h"
#include "core/include/constants.h"

namespace feasst {

class ModelTwoBody : public Model {
 public:
  double compute(
      const int group_index,
      Configuration * config,
      VisitModel * visitor) const override {
    visitor->compute(*this, config, group_index);
    return visitor->energy();
  }
  double compute(
      Configuration * config,
      VisitModel * visitor) const override {
    visitor->compute(*this, config, 0);
    return visitor->energy();
  }
  double compute(
      const Select& selection,
      const int group_index,
      Configuration * config,
      VisitModel * visitor) const override {
    visitor->compute(*this, selection, config, group_index);
    return visitor->energy();
  }
  double compute(
      const Select& selection,
      Configuration * config,
      VisitModel * visitor) const override {
    visitor->compute(*this, selection, config, 0);
    return visitor->energy();
  }
  virtual double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) const = 0;
  virtual ~ModelTwoBody() {}
};

class ModelTwoBodyFactory : public ModelTwoBody {
 public:
  void add_model(std::shared_ptr<ModelTwoBody> model) {
    models_.push_back(model); }
  const std::vector<std::shared_ptr<ModelTwoBody> >& models() const {
    return models_; }

  double energy(
      const double squared_distance,
      const int type1,
      const int type2,
      const ModelParams& model_params) const {
    double en = 0;
    for (const std::shared_ptr<ModelTwoBody> model : models_) {
      en += model->energy(squared_distance, type1, type2, model_params);
    }
//    int index = 0;
//    while ((index < static_cast<int>(models_.size())) &&
//           (en < NEAR_INFINITY)) {
//      en += models_[index]->energy(squared_distance,
//                                   type1, type2, model_params);
//      ++index;
//    }
    //INFO("en " << en);
    return en;
  }

 private:
  std::vector<std::shared_ptr<ModelTwoBody> > models_;
};

}  // namespace feasst

#endif  // FEASST_CORE_MODEL_TWO_BODY_H_
